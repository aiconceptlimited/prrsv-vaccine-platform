#!/usr/bin/env python3
"""
VAXINTAIC – Stage 04 : Current Public Epitope Ranking
Ranks epitopes based on a composite “Epitope Prioritization Index” (EPI)
that integrates hydrophilicity, flexibility, and conservation scores.

Outputs:
  • data/candidates/top_epitopes_public_current.csv
  • data/candidates/epitope_ranking_summary_public_current.csv
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np

# ────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
DATA_FILE = BASE_DIR / "data" / "epitopes" / "epitope_predictions_public_current.csv"
OUT_DIR = BASE_DIR / "data" / "candidates"
os.makedirs(OUT_DIR, exist_ok=True)

if not os.path.exists(DATA_FILE):
    raise SystemExit("❌ Missing epitope_predictions_public_current.csv — run Stage 03 first.")

# ────────────────────────────────────────────────
# LOAD DATA
# ────────────────────────────────────────────────
df = pd.read_csv(DATA_FILE)
print(f"📥 Loaded {len(df)} epitope records from {DATA_FILE}")

# make sure prrsv_type exists
if "prrsv_type" not in df.columns:
    raise SystemExit("❌ No prrsv_type column found in epitope predictions.")

# ────────────────────────────────────────────────
# NORMALIZE SCORES
# ────────────────────────────────────────────────
for col in ["hydrophilicity", "flexibility", "conservation", "epitope_score"]:
    if col in df.columns:
        min_v, max_v = df[col].min(), df[col].max()
        if max_v > min_v:  # avoid div by zero
            df[col] = (df[col] - min_v) / (max_v - min_v + 1e-8)

# ────────────────────────────────────────────────
# COMPOSITE SCORING (Epitope Prioritization Index)
# ────────────────────────────────────────────────
df["EPI"] = (
    0.45 * df["conservation"] +
    0.25 * df["hydrophilicity"] +
    0.20 * df["flexibility"] +
    0.10 * df["epitope_score"]
)

# ────────────────────────────────────────────────
# RANK + SELECT TOP-N PER TYPE
# ────────────────────────────────────────────────
TOP_N = 15
top_epitopes = (
    df.groupby("prrsv_type", group_keys=False)
      .apply(lambda g: g.nlargest(TOP_N, "EPI"))
      .reset_index(drop=True)
      .sort_values(["prrsv_type", "EPI"], ascending=[True, False])
)

# ────────────────────────────────────────────────
# SAVE RESULTS
# ────────────────────────────────────────────────
top_file = os.path.join(OUT_DIR, "top_epitopes_public_current.csv")
top_epitopes.to_csv(top_file, index=False)

summary = (
    df.groupby("prrsv_type")
      .agg(
          total_epitopes=("EPI", "count"),
          mean_EPI=("EPI", "mean"),
          top_EPI=("EPI", "max"),
          median_conservation=("conservation", "median")
      )
      .reset_index()
)
summary_file = os.path.join(OUT_DIR, "epitope_ranking_summary_public_current.csv")
summary.to_csv(summary_file, index=False)

# ────────────────────────────────────────────────
# PRINT SUMMARY
# ────────────────────────────────────────────────
print("\n📊 Epitope Ranking Summary:")
print(summary)
print(f"\n🏆 Top {TOP_N} epitopes per PRRSV type saved → {top_file}")
print("🧬 Stage 04 completed successfully.")

