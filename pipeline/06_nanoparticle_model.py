
# SCIENTIFIC VALIDITY NOTICE
#
# This is a computational/modelled nanoparticle/LNP delivery score.
# It is not experimental evidence of nanoparticle delivery or vaccine efficacy.
DELIVERY_MODEL_STATUS = "in_silico_modelled"

#!/usr/bin/env python3
"""
VAXINTAIC – Stage 06 : Current Public In-Silico Delivery Modeling

Computationally models LNP-related delivery scores for designed mRNA constructs.
The regression model is trained on synthetic data and is not experimentally validated.

Inputs:
  • data/mrna/mrna_construct_metadata_public_current.csv

Outputs:
  • data/models/nanoparticle_model_public_current.csv
  • data/models/nanoparticle_summary_public_current.csv
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor

# ────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
META_FILE = BASE_DIR / "data" / "mrna" / "mrna_construct_metadata_public_current.csv"
OUT_DIR = BASE_DIR / "data" / "models"
os.makedirs(OUT_DIR, exist_ok=True)

if not os.path.exists(META_FILE):
    raise SystemExit("❌ Missing mRNA construct metadata — run Stage 05 first.")

# ────────────────────────────────────────────────
# LOAD CONSTRUCT METADATA
# ────────────────────────────────────────────────
meta = pd.read_csv(META_FILE)
print(f"📥 Loaded {len(meta)} current-public mRNA constructs from {META_FILE}")

# ────────────────────────────────────────────────
# FEATURE SIMULATION
# ────────────────────────────────────────────────
# Add physicochemical features that affect encapsulation
meta["length_kb"] = meta["length_nt"] / 1000.0
meta["stability_index"] = np.clip((meta["GC_percent"] / 100) * meta["CAI"], 0, 1)

# Synthetic training data (simulated previous constructs)
np.random.seed(42)
train_X = np.random.rand(100, 3)
train_y = 0.6 * train_X[:, 0] + 0.3 * train_X[:, 1] + 0.1 * train_X[:, 2] + np.random.normal(0, 0.05, 100)

# ────────────────────────────────────────────────
# MODEL TRAINING (Random Forest)
# ────────────────────────────────────────────────
model = RandomForestRegressor(n_estimators=150, random_state=42)
model.fit(train_X, train_y)

# Prepare features for prediction
X_pred = meta[["GC_percent", "CAI", "length_kb"]].copy()
X_pred = (X_pred - X_pred.min()) / (X_pred.max() - X_pred.min() + 1e-8)  # normalize

# Predict encapsulation efficiency and delivery stability
meta["encapsulation_efficiency"] = np.clip(model.predict(X_pred.to_numpy()), 0.65, 0.98)
meta["delivery_stability"] = np.clip(0.8 * meta["encapsulation_efficiency"] + 0.2 * meta["stability_index"], 0, 1)

# Composite performance index
meta["delivery_index"] = np.round(
    (0.6 * meta["encapsulation_efficiency"] + 0.4 * meta["delivery_stability"]), 3
)

# ────────────────────────────────────────────────
# SAVE RESULTS
# ────────────────────────────────────────────────
out_csv = os.path.join(OUT_DIR, "nanoparticle_model_public_current.csv")
summary_csv = os.path.join(OUT_DIR, "nanoparticle_summary_public_current.csv")

meta.to_csv(out_csv, index=False)
meta_summary = meta.groupby("prrsv_type").agg(
    mean_encapsulation=("encapsulation_efficiency", "mean"),
    mean_stability=("delivery_stability", "mean"),
    mean_delivery_index=("delivery_index", "mean")
).reset_index()
meta_summary.to_csv(summary_csv, index=False)

# ────────────────────────────────────────────────
# PRINT SUMMARY
# ────────────────────────────────────────────────
print("\n📊 Nanoparticle Encapsulation Summary:")
print(meta_summary.to_string(index=False))
print(f"\n✅ Model output saved → {out_csv}")
print(f"📄 Summary → {summary_csv}")
print("🧬 Stage 06 completed successfully.")

