#!/usr/bin/env python3

from __future__ import annotations

import os
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

try:
    from run_provenance import new_run_id, write_status
except ImportError:
    write_status = None


# ============================================================================
# VAXINTAIC AUTHORITATIVE PIPELINE RUNNER
# ============================================================================
#
# This runner is intentionally limited to pipeline execution.
# Database integration is handled separately during the database remediation.
#
# Safety properties:
#   - absolute project paths
#   - explicit production Python interpreter
#   - fixed authoritative stage order
#   - per-stage timeout
#   - process-group termination on timeout
#   - unique run ID
#   - run-specific log
#   - nonzero exit on failure
# ============================================================================


BASE_DIR = Path(__file__).resolve().parents[1]
PIPELINE_DIR = BASE_DIR / "pipeline"
DATA_DIR = BASE_DIR / "data"
LOG_DIR = DATA_DIR / "logs"

PYTHON = Path(sys.executable).resolve()

# Maximum allowed runtime for an individual pipeline stage.
# This is deliberately configurable through an environment variable.
DEFAULT_STAGE_TIMEOUT_SECONDS = 30 * 60
STAGE_TIMEOUT_SECONDS = int(
    os.environ.get(
        "VAXINTAIC_STAGE_TIMEOUT_SECONDS",
        DEFAULT_STAGE_TIMEOUT_SECONDS,
    )
)

RUN_ID = new_run_id() if write_status else datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
RUN_LOG = LOG_DIR / f"pipeline_{RUN_ID}.log"

# ============================================================================
# AUTHORITATIVE V3 STAGE ORDER
# ============================================================================
#
# Dependency correction:
#   Stage 10 produces immunogenicity_scored.csv
#   Stage 07 consumes immunogenicity_scored.csv
#
# Therefore Stage 10 MUST precede Stage 07.
#
# Legacy Stage 08 and Stage 12 are intentionally excluded until their
# legacy/v3 filename and dependency relationships are reconciled separately.
# ============================================================================


PROVENANCE_ARTIFACTS = [
    "data/sequences/prrsv_orf5or6_clean.fasta",
    "data/sequences/metadata_v3.csv",
    "data/sequences/accession_manifest.csv",
    "data/epitopes/epitope_predictions_public_current.csv",
    "data/candidates/top_epitopes_public_current.csv",
    "data/mrna/mrna_constructs_public_current.fasta",
    "data/mrna/mrna_construct_metadata_public_current.csv",
    "data/models/nanoparticle_model_public_current.csv",
    "data/models/nanoparticle_summary_public_current.csv",
    "data/epitopes/immunogenicity_scored_public_current.csv",
    "data/reports/model_scores_public_current.csv",
    "data/reports/prioritization_summary_public_current.csv",
    "data/mrna/codon_adaptation_public_current.csv",
]

STAGES = [
    "01_fetch_sequences.py",
    "02_align_sequences.py",
    "03_epitope_prediction.py",
    "03b_plot_epitope_landscape.py",
    "04_epitope_ranking.py",
    "05_mrna_design.py",
    "06_nanoparticle_model.py",
    "10_immunogenicity_score.py",
    "07_in_silico_prioritization.py",
    "09_diversity_analysis.py",
    "11_codon_adaptation.py",
]


if write_status:
    write_status(RUN_ID, "STARTED", artifacts=PROVENANCE_ARTIFACTS)


def log(message: str) -> None:
    """Write a timestamped message to stdout and the run-specific log."""
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    line = f"{timestamp} [{RUN_ID}] {message}"

    print(line, flush=True)

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with RUN_LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def validate_preflight() -> None:
    """Validate runner prerequisites before executing scientific stages."""
    if not BASE_DIR.is_dir():
        raise RuntimeError(f"Missing project directory: {BASE_DIR}")

    if not PIPELINE_DIR.is_dir():
        raise RuntimeError(f"Missing pipeline directory: {PIPELINE_DIR}")

    if not DATA_DIR.is_dir():
        raise RuntimeError(f"Missing data directory: {DATA_DIR}")

    if not PYTHON.is_file():
        raise RuntimeError(f"Python interpreter not found: {PYTHON}")

    if not os.access(PYTHON, os.X_OK):
        raise RuntimeError(f"Python interpreter is not executable: {PYTHON}")

    missing = [
        stage
        for stage in STAGES
        if not (PIPELINE_DIR / stage).is_file()
    ]

    if missing:
        raise RuntimeError(
            "Missing required pipeline stages: " + ", ".join(missing)
        )

    if STAGE_TIMEOUT_SECONDS <= 0:
        raise RuntimeError(
            "VAXINTAIC_STAGE_TIMEOUT_SECONDS must be greater than zero"
        )

    log(f"Python interpreter: {PYTHON}")
    log(f"Project directory: {BASE_DIR}")
    log(f"Stage timeout: {STAGE_TIMEOUT_SECONDS} seconds")
    log(f"Authoritative stage count: {len(STAGES)}")


def terminate_process_group(process: subprocess.Popen) -> None:
    """Terminate a timed-out stage and its child processes."""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return

    try:
        process.wait(timeout=10)
        return
    except subprocess.TimeoutExpired:
        pass

    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass

    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass


def run_stage(stage_number: int, stage: str) -> float:
    """Execute one pipeline stage with timeout and process-group isolation."""
    stage_path = PIPELINE_DIR / stage

    log(
        f"START stage {stage_number:02d}/{len(STAGES):02d}: "
        f"{stage}"
    )

    start = time.monotonic()

    process = subprocess.Popen(
        [str(PYTHON), str(stage_path)],
        cwd=str(BASE_DIR),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        start_new_session=True,
    )

    try:
        stdout, _ = process.communicate(timeout=STAGE_TIMEOUT_SECONDS)

    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - start

        log(
            f"TIMEOUT stage {stage}: exceeded "
            f"{STAGE_TIMEOUT_SECONDS} seconds "
            f"(elapsed {elapsed:.1f}s)"
        )

        terminate_process_group(process)

        # Collect any output remaining after termination without imposing
        # another long wait on the timed-out stage.
        try:
            stdout, _ = process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            stdout = ""

        if stdout:
            print(stdout, end="", flush=True)
            with RUN_LOG.open("a", encoding="utf-8") as handle:
                handle.write(stdout)

        raise RuntimeError(
            f"Pipeline stage timed out: {stage} "
            f"after {STAGE_TIMEOUT_SECONDS} seconds"
        )

    if stdout:
        print(stdout, end="", flush=True)
        with RUN_LOG.open("a", encoding="utf-8") as handle:
            handle.write(stdout)

    return_code = process.returncode
    elapsed = time.monotonic() - start

    if return_code != 0:
        log(
            f"FAIL stage {stage}: exit code {return_code} "
            f"(elapsed {elapsed:.1f}s)"
        )
        raise RuntimeError(
            f"Pipeline failed at {stage} with exit code {return_code}"
        )

    log(
        f"PASS stage {stage}: "
        f"elapsed {elapsed:.1f}s"
    )

    return elapsed


def main() -> int:
    """Run the authoritative VAXINTAIC pipeline."""
    try:
        validate_preflight()

        log("=" * 72)
        log("VAXINTAIC V3 PIPELINE STARTED")
        log("=" * 72)

        total_start = time.monotonic()

        for stage_number, stage in enumerate(STAGES, 1):
            run_stage(stage_number, stage)

        total_elapsed = time.monotonic() - total_start

        log("=" * 72)
        log(
            f"VAXINTAIC V3 PIPELINE COMPLETED SUCCESSFULLY "
            f"in {total_elapsed:.1f}s"
        )
        log(f"Run log: {RUN_LOG}")
        log("=" * 72)
        if write_status:
            write_status(RUN_ID, "SUCCESS", artifacts=PROVENANCE_ARTIFACTS)
        return 0

    except KeyboardInterrupt:
        log("PIPELINE INTERRUPTED BY USER")
        if write_status:
            write_status(RUN_ID, "INTERRUPTED", error="KeyboardInterrupt", artifacts=PROVENANCE_ARTIFACTS)
        return 130

    except Exception as exc:
        log(f"PIPELINE FAILED: {exc}")
        log(f"Run log: {RUN_LOG}")
        if write_status:
            write_status(RUN_ID, "FAILED", error=str(exc), artifacts=PROVENANCE_ARTIFACTS)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
