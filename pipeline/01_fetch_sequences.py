#!/usr/bin/env python3
"""
VAXINTAIC v3 — Stage 01: Reproducible PRRSV sequence acquisition.

Default mode:
    frozen
    Uses data/sequences/accession_manifest.csv and the existing local FASTA
    as the reproducibility baseline. No network access is required.

Explicit acquisition mode:
    VAXINTAIC_ACQUISITION_MODE=ncbi
    Queries NCBI, retrieves the selected accessions, validates them, and
    writes a new accession manifest before replacing the working outputs.

The current frozen dataset is intentionally preserved unless explicit NCBI
acquisition mode is requested.
"""

from __future__ import annotations

import csv
import datetime as dt
import os
import random
import tempfile
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
SEQUENCE_DIR = DATA_DIR / "sequences"
REFERENCE_DIR = DATA_DIR / "references"

MANIFEST_FILE = SEQUENCE_DIR / "accession_manifest.csv"
FASTA_FILE = SEQUENCE_DIR / "prrsv_orf5or6_clean.fasta"
METADATA_FILE = SEQUENCE_DIR / "metadata_v3.csv"
FETCH_LOG_FILE = SEQUENCE_DIR / "fetch_log.csv"

ACQUISITION_MODE = os.environ.get(
    "VAXINTAIC_ACQUISITION_MODE", "frozen"
).strip().lower()

MAX_RECORDS = 80
MIN_LENGTH = 600
MAX_LENGTH = 750

TYPE_I_QUERY = (
    "PRRSV-1[All Fields] AND "
    "(ORF5[Gene] OR ORF6[Gene]) AND "
    "600:750[Sequence Length]"
)

TYPE_II_QUERY = (
    "PRRSV-2[All Fields] AND "
    "(ORF5[Gene] OR ORF6[Gene]) AND "
    "600:750[Sequence Length]"
)

ENTREZ_EMAIL = os.environ.get(
    "VAXINTAIC_ENTREZ_EMAIL",
    "vaxintaic@ai-lab.org",
)

MAX_RETRIES = int(os.environ.get("VAXINTAIC_NCBI_MAX_RETRIES", "5"))
BASE_BACKOFF_SECONDS = float(
    os.environ.get("VAXINTAIC_NCBI_BACKOFF_SECONDS", "5")
)


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )


def atomic_write_dataframe(df: pd.DataFrame, destination: Path) -> None:
    """Write a CSV atomically so failed writes do not corrupt the target."""
    destination.parent.mkdir(parents=True, exist_ok=True)

    fd, temp_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
        text=True,
    )

    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            df.to_csv(handle, index=False)
            handle.flush()
            os.fsync(handle.fileno())

        os.replace(temp_name, destination)

    except Exception:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def atomic_write_fasta(records, destination: Path) -> None:
    """Write FASTA atomically."""
    destination.parent.mkdir(parents=True, exist_ok=True)

    fd, temp_name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp",
        dir=destination.parent,
        text=True,
    )

    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            SeqIO.write(records, handle, "fasta")
            handle.flush()
            os.fsync(handle.fileno())

        os.replace(temp_name, destination)

    except Exception:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def validate_sequence_record(record) -> None:
    """Validate one nucleotide sequence using the standard IUPAC DNA alphabet."""
    sequence = str(record.seq).upper()

    if not sequence:
        raise ValueError(f"{record.id}: empty sequence")

    allowed_symbols = set("ACGTRYSWKMBDHVN")
    unexpected = sorted(set(sequence) - allowed_symbols)

    if unexpected:
        raise ValueError(
            f"{record.id}: sequence contains unexpected nucleotide symbols: "
            f"{unexpected}"
        )

    if not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:
        raise ValueError(
            f"{record.id}: length {len(sequence)} outside "
            f"{MIN_LENGTH}-{MAX_LENGTH} nt"
        )


def validate_records(records) -> None:
    """Validate IDs, uniqueness, sequence content, and lengths."""
    if not records:
        raise ValueError("No sequence records available")

    ids = [str(record.id) for record in records]

    if len(ids) != len(set(ids)):
        duplicates = sorted(
            {identifier for identifier in ids if ids.count(identifier) > 1}
        )
        raise ValueError(
            f"Duplicate sequence IDs detected: {duplicates[:10]}"
        )

    for record in records:
        validate_sequence_record(record)


def load_manifest() -> pd.DataFrame:
    """Load and validate the accession manifest."""
    if not MANIFEST_FILE.exists():
        raise FileNotFoundError(
            f"Accession manifest not found: {MANIFEST_FILE}"
        )

    manifest = pd.read_csv(MANIFEST_FILE, dtype=str)

    required = {
        "manifest_order",
        "id",
        "source_database",
        "acquisition_mode",
        "prrsv_type",
    }

    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(
            f"Manifest missing required columns: {sorted(missing)}"
        )

    if manifest.empty:
        raise ValueError("Accession manifest is empty")

    manifest["id"] = manifest["id"].astype(str).str.strip()
    manifest["prrsv_type"] = (
        manifest["prrsv_type"].astype(str).str.strip()
    )

    if manifest["id"].eq("").any():
        raise ValueError("Manifest contains empty accession IDs")

    if manifest["id"].duplicated().any():
        duplicates = manifest.loc[
            manifest["id"].duplicated(), "id"
        ].tolist()
        raise ValueError(
            f"Manifest contains duplicate accessions: {duplicates[:10]}"
        )

    allowed_types = {"Type I", "Type II"}

    invalid_types = sorted(
        set(manifest["prrsv_type"]) - allowed_types
    )

    if invalid_types:
        raise ValueError(
            f"Manifest contains invalid PRRSV types: {invalid_types}"
        )

    try:
        manifest["manifest_order"] = manifest["manifest_order"].astype(int)
    except ValueError as exc:
        raise ValueError(
            "manifest_order must contain integers"
        ) from exc

    if manifest["manifest_order"].duplicated().any():
        raise ValueError("Manifest contains duplicate manifest_order values")

    manifest = manifest.sort_values(
        "manifest_order",
        kind="stable",
    ).reset_index(drop=True)

    expected_order = list(range(1, len(manifest) + 1))

    if manifest["manifest_order"].tolist() != expected_order:
        raise ValueError(
            "Manifest order must be contiguous starting at 1"
        )

    return manifest


def manifest_from_local_dataset() -> pd.DataFrame:
    """
    Validate that the frozen manifest corresponds exactly to the local FASTA.
    """
    manifest = load_manifest()

    if not FASTA_FILE.exists():
        raise FileNotFoundError(
            f"Frozen FASTA not found: {FASTA_FILE}"
        )

    records = list(SeqIO.parse(FASTA_FILE, "fasta"))
    validate_records(records)

    fasta_ids = [str(record.id) for record in records]
    manifest_ids = manifest["id"].tolist()

    if set(fasta_ids) != set(manifest_ids):
        missing_from_fasta = sorted(set(manifest_ids) - set(fasta_ids))
        extra_in_fasta = sorted(set(fasta_ids) - set(manifest_ids))

        raise ValueError(
            "Frozen manifest and FASTA do not contain the same accessions. "
            f"Missing from FASTA: {missing_from_fasta[:10]}; "
            f"extra in FASTA: {extra_in_fasta[:10]}"
        )

    if len(records) != len(manifest):
        raise ValueError(
            f"Manifest has {len(manifest)} records but FASTA has "
            f"{len(records)} records"
        )

    # Reorder local records exactly according to manifest_order.
    record_by_id = {str(record.id): record for record in records}
    ordered_records = [
        record_by_id[accession]
        for accession in manifest_ids
    ]

    return manifest, ordered_records


# ---------------------------------------------------------------------------
# NCBI acquisition
# ---------------------------------------------------------------------------

def ncbi_call(function, *args, **kwargs):
    """
    Execute an Entrez request with bounded exponential backoff and jitter.
    """
    last_error = None

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return function(*args, **kwargs)

        except Exception as exc:
            last_error = exc

            if attempt == MAX_RETRIES:
                break

            delay = (
                BASE_BACKOFF_SECONDS * (2 ** (attempt - 1))
                + random.uniform(0, 1)
            )

            print(
                f"⚠️ NCBI request failed "
                f"(attempt {attempt}/{MAX_RETRIES}): {exc}"
            )
            print(f"   Retrying in {delay:.1f} seconds...")
            time.sleep(delay)

    raise RuntimeError(
        f"NCBI request failed after {MAX_RETRIES} attempts"
    ) from last_error


def ncbi_search(query: str, limit: int) -> list[str]:
    """Search NCBI and return nucleotide accessions."""
    handle = ncbi_call(
        Entrez.esearch,
        db="nucleotide",
        term=query,
        retmax=limit,
        sort="accession",
    )

    try:
        result = Entrez.read(handle)
    finally:
        handle.close()

    ids = [str(value) for value in result.get("IdList", [])]

    if not ids:
        return []

    # Convert NCBI internal IDs to accession identifiers by fetching FASTA.
    handle = ncbi_call(
        Entrez.efetch,
        db="nucleotide",
        id=ids,
        rettype="acc",
        retmode="text",
    )

    try:
        accessions = [
            line.strip()
            for line in handle.read().splitlines()
            if line.strip()
        ]
    finally:
        handle.close()

    return accessions


def ncbi_fetch_records(accessions: list[str]):
    """Fetch complete FASTA records for the requested accessions."""
    if not accessions:
        return []

    handle = ncbi_call(
        Entrez.efetch,
        db="nucleotide",
        id=accessions,
        rettype="fasta",
        retmode="text",
    )

    try:
        records = list(SeqIO.parse(handle, "fasta"))
    finally:
        handle.close()

    return records


def acquire_from_ncbi() -> tuple[pd.DataFrame, list]:
    """
    Perform an explicit NCBI acquisition and build a new manifest.
    """
    Entrez.email = ENTREZ_EMAIL

    print("🌐 NCBI acquisition mode enabled")
    print(f"   Type I query:  {TYPE_I_QUERY}")
    print(f"   Type II query: {TYPE_II_QUERY}")

    per_type = MAX_RECORDS // 2

    type_i_accessions = ncbi_search(TYPE_I_QUERY, per_type)
    type_ii_accessions = ncbi_search(TYPE_II_QUERY, per_type)

    if not type_i_accessions and not type_ii_accessions:
        raise RuntimeError("NCBI returned no accessions")

    type_i_records = ncbi_fetch_records(type_i_accessions)
    type_ii_records = ncbi_fetch_records(type_ii_accessions)

    for record in type_i_records:
        record.annotations["forced_type"] = "Type I"

    for record in type_ii_records:
        record.annotations["forced_type"] = "Type II"

    records = type_i_records + type_ii_records

    validate_records(records)

    if len(records) > MAX_RECORDS:
        raise ValueError(
            f"NCBI acquisition returned {len(records)} records; "
            f"maximum allowed is {MAX_RECORDS}"
        )

    rows = []

    for order, record in enumerate(records, start=1):
        rows.append(
            {
                "manifest_order": order,
                "id": str(record.id),
                "source_database": "NCBI nucleotide",
                "acquisition_mode": "ncbi",
                "prrsv_type": record.annotations.get(
                    "forced_type", "Unknown"
                ),
            }
        )

    manifest = pd.DataFrame(rows)

    if manifest["id"].duplicated().any():
        raise ValueError("NCBI acquisition produced duplicate accessions")

    return manifest, records


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

def build_metadata(records, manifest: pd.DataFrame) -> pd.DataFrame:
    """Build metadata from the validated sequence records."""
    type_by_id = dict(
        zip(
            manifest["id"].astype(str),
            manifest["prrsv_type"].astype(str),
        )
    )

    rows = []
    timestamp = utc_now()

    for record in records:
        accession = str(record.id)
        description = str(record.description)

        rows.append(
            {
                "id": accession,
                "organism": description,
                "country": record.annotations.get(
                    "country", "N/A"
                ),
                "collection_date": record.annotations.get(
                    "date", "N/A"
                ),
                "host": record.annotations.get(
                    "host", "N/A"
                ),
                "gene": "ORF5",
                "prrsv_type": type_by_id.get(
                    accession, "Unknown"
                ),
                "length": len(record.seq),
                "date_cached": timestamp,
                "metadata_completeness": 0.75,
            }
        )

    metadata = pd.DataFrame(rows)

    if metadata.empty:
        raise ValueError("Generated metadata is empty")

    if metadata["id"].duplicated().any():
        raise ValueError("Generated metadata contains duplicate IDs")

    if set(metadata["id"]) != set(manifest["id"]):
        raise ValueError(
            "Generated metadata IDs do not match manifest IDs"
        )

    return metadata


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def append_fetch_log(
    metadata: pd.DataFrame,
    mode: str,
) -> None:
    """Append one validated acquisition/frozen-run event to fetch_log.csv."""
    type_counts = metadata["prrsv_type"].value_counts()

    row = {
        "timestamp": utc_now(),
        "source": "local_frozen_manifest"
        if mode == "frozen"
        else "NCBI",
        "num_sequences": len(metadata),
        "avg_length": float(metadata["length"].mean()),
        "avg_completeness": float(
            metadata["metadata_completeness"].mean()
        ),
        "type_I": int(type_counts.get("Type I", 0)),
        "type_II": int(type_counts.get("Type II", 0)),
        "unknown_type": int(type_counts.get("Unknown", 0)),
    }

    FETCH_LOG_FILE.parent.mkdir(parents=True, exist_ok=True)

    exists = FETCH_LOG_FILE.exists()

    with FETCH_LOG_FILE.open(
        "a",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=row.keys(),
        )

        if not exists:
            writer.writeheader()

        writer.writerow(row)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    SEQUENCE_DIR.mkdir(parents=True, exist_ok=True)
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

    if ACQUISITION_MODE not in {"frozen", "ncbi"}:
        raise SystemExit(
            "Invalid VAXINTAIC_ACQUISITION_MODE. "
            "Use 'frozen' or 'ncbi'."
        )

    print("=" * 72)
    print("VAXINTAIC Stage 01 — Reproducible Sequence Acquisition")
    print("=" * 72)
    print(f"Base directory: {BASE_DIR}")
    print(f"Acquisition mode: {ACQUISITION_MODE}")
    print(f"Manifest: {MANIFEST_FILE}")

    if ACQUISITION_MODE == "frozen":
        print("\n🔒 FROZEN MODE")
        print("   No NCBI network request will be made.")

        manifest, records = manifest_from_local_dataset()

    else:
        print("\n🌐 NCBI MODE")
        manifest, records = acquire_from_ncbi()

        # Validate the newly acquired dataset before any output replacement.
        validate_records(records)

        # Manifest is only written after acquisition and sequence validation.
        atomic_write_dataframe(
            manifest,
            MANIFEST_FILE,
        )

    metadata = build_metadata(records, manifest)

    # Final cross-check before output.
    if len(records) != len(manifest):
        raise ValueError(
            "Final record/manifest count mismatch"
        )

    if len(metadata) != len(records):
        raise ValueError(
            "Final metadata/record count mismatch"
        )

    validate_records(records)

    # Preserve manifest order.
    record_ids = [str(record.id) for record in records]
    manifest_ids = manifest["id"].astype(str).tolist()

    if record_ids != manifest_ids:
        raise ValueError(
            "Final sequence order does not match manifest order"
        )

    # In frozen mode we deliberately do not rewrite the existing FASTA or
    # metadata. They are the frozen baseline that was validated above.
    #
    # In NCBI mode, outputs are replaced only after all validation passes.
    if ACQUISITION_MODE == "ncbi":
        atomic_write_fasta(records, FASTA_FILE)
        atomic_write_dataframe(metadata, METADATA_FILE)

    append_fetch_log(metadata, ACQUISITION_MODE)

    print("\n=== STAGE 01 VALIDATION PASSED ===")
    print(f"Records: {len(records)}")
    print(
        "Type I:",
        int((manifest["prrsv_type"] == "Type I").sum()),
    )
    print(
        "Type II:",
        int((manifest["prrsv_type"] == "Type II").sum()),
    )
    print(
        "Length range:",
        f"{metadata['length'].min()}-{metadata['length'].max()} nt",
    )
    print(
        "Unique accessions:",
        metadata["id"].nunique(),
    )

    print("\nStage 01 completed successfully.")


if __name__ == "__main__":
    main()
