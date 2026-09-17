#!/usr/bin/env python3
"""
VAXINTAIC v3 — Stage 11
Host-specific codon adaptation analysis.

FAIL-CLOSED:
CAI is not reported unless a documented Sus scrofa codon-use
reference is available.
"""

from pathlib import Path
import pandas as pd
from Bio import SeqIO
import math
import sys

BASE_DIR = Path(__file__).resolve().parents[1]
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))
INPUT_FASTA = BASE_DIR / "data" / "mrna" / "mrna_constructs_public_current.fasta"
OUTPUT_CSV = BASE_DIR / "data" / "mrna" / "codon_adaptation_public_current.csv"

from config.codon_reference import (
    TARGET_HOST_ORGANISM,
    SUS_SCROFA_CODON_USAGE,
    CODON_REFERENCE_STATUS,
    CODON_REFERENCE_SOURCE,
)

if TARGET_HOST_ORGANISM != "Sus scrofa":
    raise SystemExit(
        f"Unsupported target host for this validated configuration: "
        f"{TARGET_HOST_ORGANISM}"
    )

if not SUS_SCROFA_CODON_USAGE:
    raise SystemExit(
        "CAI BLOCKED: no documented Sus scrofa codon-usage reference is installed. "
        "Do not substitute E. coli or an invented frequency table."
    )

if not INPUT_FASTA.exists():
    raise SystemExit(f"Missing mRNA FASTA: {INPUT_FASTA}")

records = list(SeqIO.parse(INPUT_FASTA, "fasta"))

def compute_cai(sequence, ref):
    # The frozen Kazusa reference uses RNA codon notation (U).
    # Normalize the mRNA sequence to the same representation before lookup.
    sequence = sequence.upper().replace("T", "U")

    codons = [
        sequence[i:i+3]
        for i in range(0, len(sequence) - 2, 3)
    ]

    stop_codons = {"UAA", "UAG", "UGA"}

    # CAI uses the 61 sense codons. Stop codons are excluded from
    # the geometric mean because they have no relative-adaptiveness weight.
    sense_codons = [c for c in codons if c not in stop_codons]

    if not sense_codons:
        raise ValueError("sequence contains no sense codons for CAI calculation")

    weights = [ref.get(c) for c in sense_codons]

    if any(x is None or x <= 0 for x in weights):
        missing = sorted({
            c for c, x in zip(sense_codons, weights)
            if x is None or x <= 0
        })
        raise ValueError(
            "sequence contains sense codons absent from reference: "
            + ", ".join(missing)
        )

    # Geometric mean of relative codon weights.
    return round(
        math.exp(sum(math.log(x) for x in weights) / len(weights)),
        4,
    )

rows = []

for rec in records:
    seq = str(rec.seq).upper().replace("U", "T")

    if len(seq) % 3 != 0:
        raise SystemExit(f"{rec.id}: coding sequence length is not divisible by 3")

    cai = compute_cai(seq, SUS_SCROFA_CODON_USAGE)
    gc = round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2)

    rows.append({
        "construct_id": rec.id,
        "sequence_length": len(seq),
        "GC_percent": gc,
        "CAI": cai,
        "host_organism": TARGET_HOST_ORGANISM,
        "codon_reference_status": CODON_REFERENCE_STATUS,
        "codon_reference_source": CODON_REFERENCE_SOURCE,
    })

pd.DataFrame(rows).to_csv(OUTPUT_CSV, index=False)

print(f"Host: {TARGET_HOST_ORGANISM}")
print(f"Reference: {CODON_REFERENCE_SOURCE}")
print(f"CAI results saved: {OUTPUT_CSV}")
