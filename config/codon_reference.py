"""
VAXINTAIC validated Sus scrofa codon-usage reference.

Reference derived from the Kazusa codon-usage database:
  organism: Sus scrofa
  taxon: 9823
  CDS records: 2953
  sense codons: 61

The frozen JSON reference and provenance files are retained under:
  data/references/codon_usage/kazusa_9823/

The reference must not be silently replaced by another organism's codon usage.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parents[1]

REFERENCE_PATH = (
    BASE_DIR
    / "data"
    / "references"
    / "codon_usage"
    / "kazusa_9823"
    / "sus_scrofa_9823_cai_reference.json"
)

PROVENANCE_PATH = (
    BASE_DIR
    / "data"
    / "references"
    / "codon_usage"
    / "kazusa_9823"
    / "provenance.json"
)

TARGET_HOST_ORGANISM = "Sus scrofa"
TARGET_TAXON_ID = "9823"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


if not REFERENCE_PATH.is_file():
    raise RuntimeError(
        f"Missing validated codon reference: {REFERENCE_PATH}"
    )

if not PROVENANCE_PATH.is_file():
    raise RuntimeError(
        f"Missing codon-reference provenance: {PROVENANCE_PATH}"
    )


REFERENCE_SHA256 = _sha256(REFERENCE_PATH)
PROVENANCE_SHA256 = _sha256(PROVENANCE_PATH)

with REFERENCE_PATH.open("r", encoding="utf-8") as fh:
    _reference = json.load(fh)

with PROVENANCE_PATH.open("r", encoding="utf-8") as fh:
    PROVENANCE = json.load(fh)


if PROVENANCE.get("organism") != TARGET_HOST_ORGANISM:
    raise RuntimeError(
        "SCIENTIFIC VALIDITY ERROR: codon reference organism mismatch."
    )

if str(PROVENANCE.get("taxon_id")) != TARGET_TAXON_ID:
    raise RuntimeError(
        "SCIENTIFIC VALIDITY ERROR: codon reference taxon mismatch."
    )

CODON_REFERENCE_SOURCE = PROVENANCE.get("source_url")
CODON_REFERENCE_STATUS = "validated_frozen_kazusa_9823"


# Accept the frozen reference whether the JSON stores the frequencies
# directly or under the explicit codon-usage mapping key.
if isinstance(_reference, dict):
    if "codon_usage" in _reference:
        _raw_usage = _reference["codon_usage"]
    elif "relative_adaptiveness" in _reference:
        _raw_usage = _reference["relative_adaptiveness"]
    else:
        _raw_usage = _reference

    # The frozen Kazusa-derived JSON stores metadata for each codon:
    # {"amino_acid": "...", "count": ..., "relative_adaptiveness": ...}
    # Stage 11 needs the numerical relative-adaptiveness weights.
    SUS_SCROFA_CODON_USAGE = {
        codon: (
            record["relative_adaptiveness"]
            if isinstance(record, dict)
            else record
        )
        for codon, record in _raw_usage.items()
    }
else:
    raise RuntimeError(
        "SCIENTIFIC VALIDITY ERROR: invalid codon-reference JSON structure."
    )


if not isinstance(SUS_SCROFA_CODON_USAGE, dict):
    raise RuntimeError(
        "SCIENTIFIC VALIDITY ERROR: codon reference must be a mapping."
    )

if len(SUS_SCROFA_CODON_USAGE) != 61:
    raise RuntimeError(
        "SCIENTIFIC VALIDITY ERROR: expected exactly 61 sense-codon weights."
    )

STOP_CODONS = {"UAA", "UAG", "UGA"}

if STOP_CODONS.intersection(SUS_SCROFA_CODON_USAGE):
    raise RuntimeError(
        "SCIENTIFIC VALIDITY ERROR: stop codon found in CAI reference."
    )

if set(SUS_SCROFA_CODON_USAGE) & STOP_CODONS:
    raise RuntimeError(
        "SCIENTIFIC VALIDITY ERROR: stop codons cannot have CAI weights."
    )

for codon, weight in SUS_SCROFA_CODON_USAGE.items():
    if not isinstance(codon, str) or len(codon) != 3:
        raise RuntimeError(
            f"SCIENTIFIC VALIDITY ERROR: invalid codon label: {codon!r}"
        )
    if not isinstance(weight, (int, float)) or not (0 < float(weight) <= 1):
        raise RuntimeError(
            f"SCIENTIFIC VALIDITY ERROR: invalid CAI weight for {codon}: {weight!r}"
        )

# Build the amino-acid -> synonymous-codon mapping required by mRNA design.
# The authoritative numerical reference remains codon -> relative-adaptiveness.
from Bio.Data import CodonTable

_table = CodonTable.unambiguous_rna_by_name["Standard"]
SUS_SCROFA_CODONS_BY_AA = {}

for codon, amino_acid in _table.forward_table.items():
    if codon in SUS_SCROFA_CODON_USAGE:
        SUS_SCROFA_CODONS_BY_AA.setdefault(amino_acid, []).append(codon)

EXPECTED_AA = set("ACDEFGHIKLMNPQRSTVWY")

if set(SUS_SCROFA_CODONS_BY_AA) != EXPECTED_AA:
    missing = sorted(EXPECTED_AA - set(SUS_SCROFA_CODONS_BY_AA))
    extra = sorted(set(SUS_SCROFA_CODONS_BY_AA) - EXPECTED_AA)
    raise RuntimeError(
        "SCIENTIFIC VALIDITY ERROR: incomplete amino-acid codon mapping; "
        f"missing={missing}, extra={extra}"
    )

for amino_acid in SUS_SCROFA_CODONS_BY_AA:
    SUS_SCROFA_CODONS_BY_AA[amino_acid] = sorted(
        SUS_SCROFA_CODONS_BY_AA[amino_acid]
    )
