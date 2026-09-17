#!/usr/bin/env python3
"""
VAXINTAIC v3 – Stage 03b : Epitope Landscape Visualization

Reads the epitope_predictions_v3.csv file and plots
epitope score vs position for each PRRSV type.
Generates publication-quality visualizations.
"""

import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# ────────────────────────────────────────────────
# Paths
# ────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
DATA_FILE = BASE_DIR / "data" / "epitopes" / "epitope_predictions_v3.csv"
PLOT_DIR = BASE_DIR / "data" / "epitopes" / "plots"
os.makedirs(PLOT_DIR, exist_ok=True)

# ────────────────────────────────────────────────
# Load data
# ────────────────────────────────────────────────
if not os.path.exists(DATA_FILE):
    raise SystemExit("❌ Missing data/epitopes/epitope_predictions_v3.csv — run Stage 03 first.")

df = pd.read_csv(DATA_FILE)
types = df["prrsv_type"].unique().tolist()
print(f"📊 Loaded {len(df)} epitope records from: {DATA_FILE}")

# ────────────────────────────────────────────────
# Plot configuration
# ────────────────────────────────────────────────
plt.rcParams.update({
    "font.size": 11,
    "figure.figsize": (10, 5),
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True
})

# ────────────────────────────────────────────────
# Generate plots per type
# ────────────────────────────────────────────────
for t in types:
    subset = df[df["prrsv_type"] == t].copy()
    subset["midpoint"] = (subset["start"] + subset["end"]) / 2

    plt.figure()
    plt.plot(subset["midpoint"], subset["epitope_score"], linewidth=2)
    plt.title(f"Epitope Landscape – {t}")
    plt.xlabel("Position (AA)")
    plt.ylabel("Epitope Score")
    plt.ylim(0, 1)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()

    out_path = os.path.join(PLOT_DIR, f"{t.replace(' ', '_').lower()}_epitope_landscape.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"✅ Saved plot → {out_path}")

# ────────────────────────────────────────────────
# Combined overlay plot
# ────────────────────────────────────────────────
plt.figure(figsize=(10, 5))
for t in types:
    subset = df[df["prrsv_type"] == t].copy()
    subset["midpoint"] = (subset["start"] + subset["end"]) / 2
    plt.plot(subset["midpoint"], subset["epitope_score"], linewidth=2, label=t)

plt.title("PRRSV Epitope Landscape Comparison (Type I vs Type II)")
plt.xlabel("Position (AA)")
plt.ylabel("Epitope Score")
plt.ylim(0, 1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

combined_path = os.path.join(PLOT_DIR, "combined_epitope_landscape.png")
plt.savefig(combined_path, dpi=300)
plt.close()

print(f"📈 Combined comparison plot saved → {combined_path}")
print("🧬 Stage 03b completed successfully.")

