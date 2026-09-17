"""Run-level provenance and status utilities."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parents[1]
RUN_DIR = BASE_DIR / "data" / "run"
RUN_DIR.mkdir(parents=True, exist_ok=True)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def git_revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=BASE_DIR,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def new_run_id() -> str:
    return (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + "_"
        + uuid.uuid4().hex[:12]
    )


def write_status(
    run_id: str,
    status: str,
    *,
    stage: str | None = None,
    error: str | None = None,
    artifacts: list[str] | None = None,
) -> Path:
    payload = {
        "run_id": run_id,
        "status": status,
        "stage": stage,
        "error": error,
        "timestamp_utc": utc_now(),
        "git_revision": git_revision(),
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "pid": __import__("os").getpid(),
        "artifacts": {},
    }

    if artifacts:
        for relative in artifacts:
            path = BASE_DIR / relative
            if path.exists() and path.is_file():
                payload["artifacts"][relative] = {
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }

    target = RUN_DIR / f"{run_id}.json"
    tmp = target.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(target)
    return target
