#!/usr/bin/env python3
"""
VAXINTAIC v3 — Stage 03
Protein-level epitope landscape generation.

The authoritative PRRSV ORF nucleotide sequences are translated in
the established frame (frame 0 of the coding sequence). Candidate
windows are therefore amino-acid sequences, not arbitrary nucleotide
windows.

IMPORTANT:
This remains a sequence-derived computational scoring model. It is
not an MHC-binding predictor and does not claim experimental
immunogenicity.
"""

from pathlib import Path
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_FILE = BASE_DIR / "data" / "sequences" / "prrsv_orf5or6_clean.fasta"
OUT_DIR = BASE_DIR / "data" / "epitopes"
OUTPUT_FILE = OUT_DIR / "epitope_predictions_public_current.csv"

EPITOPE_PREDICTOR_STATUS = "sequence_derived_protein_window_score"
WINDOW_AA = int(os.environ.get("VAXINTAIC_EPITOPE_WINDOW_AA", "15"))

OUT_DIR.mkdir(parents=True, exist_ok=True)

if not INPUT_FILE.exists():
    raise SystemExit(f"Missing authoritative FASTA: {INPUT_FILE}")

records = list(SeqIO.parse(INPUT_FILE, "fasta"))
MANIFEST_FILE = BASE_DIR / "data" / "sequences" / "accession_manifest.csv"
manifest = pd.read_csv(MANIFEST_FILE)
required_manifest_columns = {"id", "prrsv_type"}
if not required_manifest_columns.issubset(manifest.columns):
    raise SystemExit("accession_manifest.csv is missing required columns: id, prrsv_type")
manifest_type_by_id = dict(zip(manifest["id"].astype(str), manifest["prrsv_type"].astype(str)))
if set(manifest_type_by_id) != {rec.id for rec in records}:
    raise SystemExit("accession_manifest.csv IDs do not exactly match the authoritative FASTA IDs")


if not records:
    raise SystemExit("No PRRSV sequences found.")

rows = []

def hydrophilicity_score(peptide):
    # Transparent sequence-derived heuristic.
    hydrophilic = set("DEKRHQNST")
    return round(sum(a in hydrophilic for a in peptide) / len(peptide), 4)

def flexibility_score(peptide):
    # Transparent sequence-derived heuristic.
    flexible = set("GSA")
    return round(sum(a in flexible for a in peptide) / len(peptide), 4)

for rec in records:
    dna = str(rec.seq).upper().replace("U", "T")

    if len(dna) % 3 != 0:
        raise SystemExit(
            f"Sequence {rec.id} is not divisible by 3; refusing translation."
        )

    protein = str(Seq(dna).translate())

    # Remove terminal stop only.
    if protein.endswith("*"):
        protein = protein[:-1]

    if "*" in protein:
        raise SystemExit(
            f"Sequence {rec.id} contains an internal stop codon; refusing analysis."
        )

    prrsv_type = manifest_type_by_id.get(rec.id)
    if prrsv_type not in {"Type I", "Type II"}:
        raise SystemExit(
            f"Sequence {rec.id} has no valid PRRSV type in accession_manifest.csv"
        )

    # Non-overlapping? No: use sliding windows to preserve landscape.
    for start0 in range(0, len(protein) - WINDOW_AA + 1):
        peptide = protein[start0:start0 + WINDOW_AA]

        hydro = hydrophilicity_score(peptide)
        flex = flexibility_score(peptide)

        # Conservation is calculated across the cohort below.
        rows.append({
            "prrsv_type": prrsv_type,
            "start": start0 + 1,
            "end": start0 + WINDOW_AA,
            "consensus_seq": peptide,
            "sequence_type": "amino_acid",
            "hydrophilicity": hydro,
            "flexibility": flex,
            "source_accession": rec.id,
        })

df = pd.DataFrame(rows)

# Cohort conservation at each type/start position.
# Defined as the frequency of the dominant valid amino-acid
# window among sequences containing that window.
# This measures sequence-derived peptide identity conservation,
# not phylogenetic conservation or evolutionary rate.
valid_for_conservation = df[
    df["consensus_seq"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")
].copy()

dominant_counts = (
    valid_for_conservation
    .groupby(["prrsv_type", "start", "end"])["consensus_seq"]
    .value_counts()
    .rename("peptide_count")
    .reset_index()
)

dominant = (
    dominant_counts
    .groupby(["prrsv_type", "start", "end"])["peptide_count"]
    .max()
    .reset_index(name="dominant_peptide_count")
)

valid_counts = (
    valid_for_conservation
    .groupby(["prrsv_type", "start", "end"])["consensus_seq"]
    .size()
    .reset_index(name="valid_sequences")
)

conservation = dominant.merge(
    valid_counts,
    on=["prrsv_type", "start", "end"],
    how="left",
)

conservation["conservation"] = (
    conservation["dominant_peptide_count"]
    / conservation["valid_sequences"]
).round(4)

df = df.merge(
    conservation[
        ["prrsv_type", "start", "end", "conservation"]
    ],
    on=["prrsv_type", "start", "end"],
    how="left",
)

df["epitope_score"] = (
    0.35 * df["hydrophilicity"]
    + 0.25 * df["flexibility"]
    + 0.40 * df["conservation"]
).round(4)

# Exclude translated windows containing unknown amino-acid residues.
# Unknown residues are not assigned an invented identity.
df = df[df["consensus_seq"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")].copy()

# Validated peptide representation generated from the established frame-0 translation.
df["sequence"] = df["consensus_seq"]
assert (df["sequence_type"] == "amino_acid").all(), "SCIENTIFIC VALIDITY ERROR: non-peptide sequence detected"
assert df["sequence"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+").all(), "SCIENTIFIC VALIDITY ERROR: invalid amino-acid sequence detected"

df = df[
    [
        "prrsv_type",
        "start",
        "end",
        "consensus_seq",
        "sequence",
        "sequence_type",
        "hydrophilicity",
        "flexibility",
        "conservation",
        "epitope_score",
    ]
]

df.to_csv(OUTPUT_FILE, index=False)

print(f"EPITOPE_PREDICTOR_STATUS={EPITOPE_PREDICTOR_STATUS}")
print(f"WINDOW_AA={WINDOW_AA}")
print(f"records={len(df)}")
print(f"output={OUTPUT_FILE}")
print("Stage 03 completed.")
