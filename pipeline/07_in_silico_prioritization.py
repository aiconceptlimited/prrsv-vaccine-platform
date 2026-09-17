# SCIENTIFIC INTERPRETATION
MODEL_STATUS = "computational_unvalidated_same_dataset"
OUTPUT_INTERPRETATION = "in_silico_prioritization_score"


# SCIENTIFIC VALIDITY NOTICE
#
# This module produces an in-silico prioritization/model score.
# Unless independently validated against an external held-out biological
# dataset, the output must NOT be interpreted as experimentally demonstrated
# vaccine efficacy.
MODEL_VALIDATION_STATUS = "computational_unvalidated_same_dataset"
OUTPUT_INTERPRETATION = "in_silico_prioritization_score"

#!/usr/bin/env python3
"""
VAXINTAIC – Stage 07 : Current Public In-Silico Prioritization

Combines molecular, immunogenicity-proxy, and delivery-model metrics
to generate an exploratory in-silico prioritization score for each construct.

Inputs:
  • data/models/nanoparticle_model_public_current.csv
  • data/epitopes/immunogenicity_scored_public_current.csv

Outputs:
  • data/reports/model_scores_public_current.csv
  • data/reports/prioritization_summary_public_current.csv
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.ensemble import GradientBoostingRegressor

# ────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
LNP_FILE = BASE_DIR / "data" / "models" / "nanoparticle_model_public_current.csv"
IMMUNO_FILE = BASE_DIR / "data" / "epitopes" / "immunogenicity_scored.csv"
OUT_DIR = BASE_DIR / "data" / "reports"
os.makedirs(OUT_DIR, exist_ok=True)

if not os.path.exists(LNP_FILE):
    raise SystemExit("❌ Missing nanoparticle_model_public_current.csv — run current-public Stage 06 first.")
if not os.path.exists(IMMUNO_FILE):
    raise SystemExit("❌ Missing immunogenicity_scored.csv — ensure Stage 10 runs before this.")

# ────────────────────────────────────────────────
# LOAD DATA
# ────────────────────────────────────────────────
lnp = pd.read_csv(LNP_FILE)
immuno = pd.read_csv(IMMUNO_FILE)

# 🧠 Add PRRSV type column if missing
if "prrsv_type" not in immuno.columns:
    if "id" in immuno.columns:
        immuno["prrsv_type"] = immuno["id"].apply(
            lambda x: "Type I" if "type_i" in str(x).lower() or "_i_" in str(x).lower()
            else ("Type II" if "type_ii" in str(x).lower() or "_ii_" in str(x).lower()
                  else "Unknown")
        )
    else:
        immuno["prrsv_type"] = "Unknown"
    print("🧩 Added inferred 'prrsv_type' column to immunogenicity data.")

# ────────────────────────────────────────────────
# COMPUTE IMMUNOGENICITY STATS
# ────────────────────────────────────────────────
immuno_summary = immuno.groupby("prrsv_type")["immunogenicity"].mean().reset_index()
immuno_summary.rename(columns={"immunogenicity": "mean_immunogenicity"}, inplace=True)

df = pd.merge(lnp, immuno_summary, on="prrsv_type", how="left")
print(f"📥 Loaded data for {len(df)} constructs")

# ────────────────────────────────────────────────
# FEATURE PREPARATION
# ────────────────────────────────────────────────
X = df[[
    "GC_percent",
    "encapsulation_efficiency",
    "delivery_stability",
    "mean_immunogenicity",
]].copy()

# Missing values are filled only for modelled inputs. CAI is not used by
# this stage; host-specific codon adaptation is handled in Stage 11.
fill_defaults = {
    "GC_percent": 50.0,
    "encapsulation_efficiency": X["encapsulation_efficiency"].mean(skipna=True),
    "delivery_stability": X["delivery_stability"].mean(skipna=True),
    "mean_immunogenicity": 0.8,
}
X = X.fillna(fill_defaults)

# Synthetic training target.
# This target is generated from the same computational inputs used for
# fitting and prediction; it is not an experimentally measured endpoint.
# The resulting score is therefore an exploratory in-silico prioritization
# score, not an experimentally validated vaccine-efficacy prediction.
np.random.seed(42)
y = (
    0.25 * X["GC_percent"] / 100
    + 0.35 * X["encapsulation_efficiency"]
    + 0.25 * X["delivery_stability"]
    + 0.15 * X["mean_immunogenicity"]
    + np.random.normal(0, 0.02, len(X))
)

# ────────────────────────────────────────────────
# MODEL TRAINING & PREDICTION
# ────────────────────────────────────────────────
model = GradientBoostingRegressor(n_estimators=150, learning_rate=0.05, random_state=42)
model.fit(X, y)
df["in_silico_model_score"] = np.clip(model.predict(X), 0, 1)

# Composite vaccine performance index
if "delivery_index" in df.columns:
    df["prioritization_index"] = np.round(
        0.5 * df["in_silico_model_score"] + 0.3 * df["delivery_index"] + 0.2 * df["mean_immunogenicity"], 3
    )
else:
    df["prioritization_index"] = np.round(
        0.7 * df["in_silico_model_score"] + 0.3 * df["mean_immunogenicity"], 3
    )

# 🧩 Ensure no NaN values remain
df["prioritization_index"] = df["prioritization_index"].fillna(df["in_silico_model_score"].round(3))

# ────────────────────────────────────────────────
# OUTPUTS
# ────────────────────────────────────────────────
eff_csv = os.path.join(OUT_DIR, "model_scores_public_current.csv")
summary_csv = os.path.join(OUT_DIR, "prioritization_summary_public_current.csv")

# Do not propagate stale/unvalidated codon-adaptation values into the
# Stage 07 artifact. CAI is unavailable until a documented Sus scrofa
# codon-usage reference is installed.
df = df.drop(columns=["CAI"], errors="ignore")

df.to_csv(eff_csv, index=False)

summary = df.groupby("prrsv_type").agg(
    mean_in_silico_model_score=("in_silico_model_score", "mean"),
    mean_prioritization_index=("prioritization_index", "mean"),
    mean_GC_percent=("GC_percent", "mean")
).reset_index()

summary.to_csv(summary_csv, index=False)

# ────────────────────────────────────────────────
# REPORT SUMMARY
# ────────────────────────────────────────────────
print("\n📊 In-Silico Prioritization Summary:")
print(summary.to_string(index=False))

print(f"\n✅ Detailed results saved → {eff_csv}")
print(f"📄 Summary saved → {summary_csv}")
print("🧬 Stage 07 current-public completed successfully.")

