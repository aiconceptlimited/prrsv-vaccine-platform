
# SCIENTIFIC VALIDITY NOTICE
#
# The legacy diversity calculation is a dominance-complement metric:
#     1 - dominant_frequency
# It must not be described as Shannon entropy, Simpson diversity, or another
# standardized diversity index unless separately implemented and validated.
DIVERSITY_METRIC_NAME = "dominance_complement"
DIVERSITY_METRIC_FORMULA = "1 - dominant_frequency"

#!/usr/bin/env python3
"""
VAXINTAIC – Stage 09 : Current Public Diversity Analysis

Analyzes sequence conservation and variability within PRRSV Type I and Type II alignments.
Produces diversity profiles and identifies conserved regions for vaccine targeting.
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import AlignIO

# ────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
ALIGN_DIR = BASE_DIR / "data" / "alignments"
OUT_DIR = BASE_DIR / "data" / "analysis"
os.makedirs(OUT_DIR, exist_ok=True)

types = {
    "Type I": os.path.join(ALIGN_DIR, "prrsv_typei_aligned.fasta"),
    "Type II": os.path.join(ALIGN_DIR, "prrsv_typeii_aligned.fasta"),
}

# ────────────────────────────────────────────────
# FUNCTION
# ────────────────────────────────────────────────
def calculate_diversity(fasta_file):
    """Calculates per-position variability index."""
    alignment = AlignIO.read(fasta_file, "fasta")
    arr = np.array([list(rec.seq) for rec in alignment])
    diversity = []
    for i in range(arr.shape[1]):
        column = arr[:, i]
        bases, counts = np.unique(column, return_counts=True)
        diversity.append(1 - (np.max(counts) / np.sum(counts)))
    return diversity

# ────────────────────────────────────────────────
# MAIN LOOP
# ────────────────────────────────────────────────
summary_rows = []

for prrsv_type, path in types.items():
    if not os.path.exists(path):
        print(f"⚠️ Missing alignment for {prrsv_type}, skipping.")
        continue

    print(f"🔍 Analyzing diversity for {prrsv_type}...")
    diversity = calculate_diversity(path)
    df = pd.DataFrame({
        "position": range(1, len(diversity) + 1),
        "diversity": diversity,
        "prrsv_type": prrsv_type
    })

    out_csv = os.path.join(OUT_DIR, f"diversity_profile_{prrsv_type.replace(' ', '').lower()}_public_current.csv")
    df.to_csv(out_csv, index=False)

    mean_div = df["diversity"].mean()
    conserved = (df["diversity"] < 0.1).sum()
    variable = (df["diversity"] > 0.3).sum()

    summary_rows.append({
        "prrsv_type": prrsv_type,
        "mean_diversity": round(mean_div, 4),
        "conserved_positions": conserved,
        "variable_positions": variable,
    })

# ────────────────────────────────────────────────
# SUMMARY OUTPUT
# ────────────────────────────────────────────────
summary_df = pd.DataFrame(summary_rows)
summary_csv = os.path.join(OUT_DIR, "diversity_summary_public_current.csv")
summary_df.to_csv(summary_csv, index=False)

print("\n📊 Diversity Analysis Summary:")
print(summary_df.to_string(index=False))
print(f"\n✅ Detailed diversity profiles saved to → {OUT_DIR}")
print("🧬 Stage 09 completed successfully.")

