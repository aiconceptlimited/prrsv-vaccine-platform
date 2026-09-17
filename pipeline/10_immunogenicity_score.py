
def _require_real_epitope_sequence(df):
    """
    Scientific validity guard:
    immunogenicity scoring must operate on the actual upstream epitope
    sequence. It must never manufacture a sequence from a numerical score.
    """
    candidates = [
        "sequence",
        "epitope_sequence",
        "epitope",
        "peptide",
        "peptide_sequence",
    ]

    for col in candidates:
        if col in df.columns:
            values = df[col].fillna("").astype(str).str.strip()
            if (values != "").all():
                return col

    raise ValueError(
        "SCIENTIFIC VALIDITY ERROR: no authoritative epitope sequence "
        "column is available. Refusing pseudo-sequence reconstruction."
    )

#!/usr/bin/env python3
"""
VAXINTAIC current-public – Stage 10 : Immunogenicity Scoring (Final Robust Version)
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np

# ────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_FILE = BASE_DIR / "data" / "candidates" / "top_epitopes_public_current.csv"
OUTPUT_FILE = BASE_DIR / "data" / "epitopes" / "immunogenicity_scored_public_current.csv"
os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

if not os.path.exists(INPUT_FILE):
    raise SystemExit("❌ Missing top_epitopes_public_current.csv — run Stage 04 first.")

# ────────────────────────────────────────────────
# LOAD EPITOPE DATA
# ────────────────────────────────────────────────
epitopes = pd.read_csv(INPUT_FILE)
print(f"📥 Loaded {len(epitopes)} epitopes for immunogenicity scoring...")

# ────────────────────────────────────────────────
REQUIRE_AUTHORITATIVE_SEQUENCE = True
sequence_column = _require_real_epitope_sequence(epitopes)
if sequence_column != "sequence":
    epitopes["sequence"] = epitopes[sequence_column].astype(str)
print(f"Using authoritative epitope sequence column: {sequence_column}")

# AMINO ACID PROPERTIES TABLE
# ────────────────────────────────────────────────
AA_PROPERTIES = {
    "A": {"hydro": 1.8, "charge": 0, "polarity": 8.1},
    "R": {"hydro": -4.5, "charge": +1, "polarity": 10.5},
    "N": {"hydro": -3.5, "charge": 0, "polarity": 11.6},
    "D": {"hydro": -3.5, "charge": -1, "polarity": 13.0},
    "C": {"hydro": 2.5, "charge": 0, "polarity": 5.5},
    "Q": {"hydro": -3.5, "charge": 0, "polarity": 10.5},
    "E": {"hydro": -3.5, "charge": -1, "polarity": 12.3},
    "G": {"hydro": -0.4, "charge": 0, "polarity": 9.0},
    "H": {"hydro": -3.2, "charge": +1, "polarity": 10.4},
    "I": {"hydro": 4.5, "charge": 0, "polarity": 5.2},
    "L": {"hydro": 3.8, "charge": 0, "polarity": 4.9},
    "K": {"hydro": -3.9, "charge": +1, "polarity": 11.3},
    "M": {"hydro": 1.9, "charge": 0, "polarity": 5.7},
    "F": {"hydro": 2.8, "charge": 0, "polarity": 5.2},
    "P": {"hydro": -1.6, "charge": 0, "polarity": 8.0},
    "S": {"hydro": -0.8, "charge": 0, "polarity": 9.2},
    "T": {"hydro": -0.7, "charge": 0, "polarity": 8.6},
    "W": {"hydro": -0.9, "charge": 0, "polarity": 5.4},
    "Y": {"hydro": -1.3, "charge": 0, "polarity": 6.2},
    "V": {"hydro": 4.2, "charge": 0, "polarity": 5.9},
}

# ────────────────────────────────────────────────
# IMMUNOGENICITY FUNCTION
# ────────────────────────────────────────────────
def compute_immunogenicity(seq):
    seq = seq.upper()
    valid_aa = [AA_PROPERTIES[a] for a in seq if a in AA_PROPERTIES]
    if not valid_aa:
        return np.nan

    hydro = np.mean([a["hydro"] for a in valid_aa])
    charge = np.mean([a["charge"] for a in valid_aa])
    polarity = np.mean([a["polarity"] for a in valid_aa])

    score = (
        0.4 * (1 - abs(hydro) / 4.5)
        + 0.3 * (1 - abs(charge) / 2.0)
        + 0.3 * (1 - abs(polarity - 8) / 8)
    )
    return np.clip(score, 0, 1)

# ────────────────────────────────────────────────
# APPLY IMMUNOGENICITY CALCULATION
# ────────────────────────────────────────────────
epitopes["immunogenicity"] = epitopes["sequence"].apply(compute_immunogenicity)
epitopes["immunogenicity"] = epitopes["immunogenicity"].fillna(0.5)

# Average by PRRSV type
if "prrsv_type" in epitopes.columns:
    summary = epitopes.groupby("prrsv_type")["immunogenicity"].mean().reset_index()
else:
    summary = pd.DataFrame({
        "prrsv_type": ["Unknown"],
        "immunogenicity": [epitopes["immunogenicity"].mean()]
    })

# ────────────────────────────────────────────────
# SAVE RESULTS
# ────────────────────────────────────────────────
epitopes.to_csv(OUTPUT_FILE, index=False)

print("\n📊 Immunogenicity Summary:")
print(summary.to_string(index=False))
print(f"\n✅ Scored epitopes saved → {OUTPUT_FILE}")
print("🧬 Stage 10 completed successfully.")

