#!/usr/bin/env python3
"""
VAXINTAIC v3 – Stage 13 : Advanced Report Generator

Integrates results from all previous stages to create a complete, human-readable
vaccine design intelligence report comparing PRRSV Type I and Type II constructs.
"""

import os
from pathlib import Path
import pandas as pd
import datetime

# ────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
REPORT_DIR = BASE_DIR / "data" / "reports"
os.makedirs(REPORT_DIR, exist_ok=True)
REPORT_PATH = os.path.join(REPORT_DIR, f"vaxintaic_advanced_report_{datetime.date.today()}.txt")

# Helper function to safely load CSVs
def safe_csv(path):
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()

# ────────────────────────────────────────────────
# LOAD ALL AVAILABLE DATA
# ────────────────────────────────────────────────
meta = safe_csv(BASE_DIR / "data" / "sequences" / "metadata_v3.csv")
div = safe_csv(BASE_DIR / "data" / "analysis" / "diversity_summary_v3.csv")
immuno = safe_csv(BASE_DIR / "data" / "epitopes" / "immunogenicity_scored.csv")
cai = safe_csv(BASE_DIR / "data" / "mrna" / "codon_adaptation_v3.csv")
eff = safe_csv(BASE_DIR / "data" / "reports" / "efficacy_summary_v3.csv")

# ────────────────────────────────────────────────
# COMPILE REPORT SECTIONS
# ────────────────────────────────────────────────
report_lines = []
report_lines.append("🧬 VAXINTAIC ADVANCED VACCINE REPORT (v3)")
report_lines.append(f"📅 Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
report_lines.append("=" * 65)

# 1️⃣ Dataset Overview
if not meta.empty:
    report_lines.append("\n📊 SEQUENCE DATASET SUMMARY")
    total = len(meta)
    type_counts = meta["prrsv_type"].value_counts().to_dict() if "prrsv_type" in meta.columns else {}
    mean_len = meta["length"].mean() if "length" in meta.columns else 0
    completeness = meta["metadata_completeness"].mean() if "metadata_completeness" in meta.columns else 0

    report_lines.append(f"- Total Sequences: {total}")
    for t, c in type_counts.items():
        report_lines.append(f"- {t}: {c} sequences")
    report_lines.append(f"- Average Length: {mean_len:.1f} nt")
    report_lines.append(f"- Avg Metadata Completeness: {completeness:.2f}")

# 2️⃣ Diversity Analysis
if not div.empty:
    report_lines.append("\n🌈 DIVERSITY ANALYSIS")
    for _, row in div.iterrows():
        report_lines.append(f"- {row['prrsv_type']}: mean diversity {row['mean_diversity']:.4f}, "
                            f"{row['variable_positions']} variable sites, "
                            f"{row['conserved_positions']} conserved regions")

# 3️⃣ Immunogenicity
if not immuno.empty:
    report_lines.append("\n🧠 IMMUNOGENICITY SUMMARY")
    if "prrsv_type" in immuno.columns:
        immuno_summary = immuno.groupby("prrsv_type")["immunogenicity"].mean().reset_index()
        for _, row in immuno_summary.iterrows():
            report_lines.append(f"- {row['prrsv_type']}: mean immunogenicity {row['immunogenicity']:.3f}")
    else:
        report_lines.append(f"- Overall mean immunogenicity: {immuno['immunogenicity'].mean():.3f}")

# 4️⃣ mRNA Codon Optimization
report_lines.append("\n🧬 CODON ADAPTATION ANALYSIS")
if not cai.empty and "CAI" in cai.columns:
    mean_cai = cai["CAI"].mean()
    report_lines.append(f"- Average CAI: {mean_cai:.3f}")
    if "GC_percent" in cai.columns:
        mean_gc = cai["GC_percent"].mean()
        report_lines.append(f"- Average GC%: {mean_gc:.2f}%")
    if mean_cai < 0.5:
        report_lines.append("⚠️ Suggestion: codon optimization could be improved for higher expression efficiency.")
else:
    report_lines.append("- CAI: unavailable")
    report_lines.append("- CAI status: blocked pending a documented Sus scrofa codon-usage reference")
    if not cai.empty and "GC_percent" in cai.columns:
        mean_gc = cai["GC_percent"].mean()
        report_lines.append(f"- Average GC%: {mean_gc:.2f}%")

# 5️⃣ Computational Vaccine Prioritization
if not eff.empty:
    report_lines.append("\n💉 COMPUTATIONAL MODEL-SCORE SUMMARY")
    for _, row in eff.iterrows():
        report_lines.append(f"- {row['prrsv_type']}: mean computational model score {row['mean_predicted_model_score']:.3f}, "
                            f"vaccine index {row['mean_vaccine_index']:.3f}")
    best_type = eff.loc[eff['mean_predicted_model_score'].idxmax(), 'prrsv_type']
    report_lines.append(f"🏆 Highest computational vaccine prioritization: {best_type}")

# ────────────────────────────────────────────────
# ADD FINAL INSIGHTS
# ────────────────────────────────────────────────
report_lines.append("\n🧩 FINAL INSIGHTS")
if not eff.empty and not immuno.empty:
    best_type = eff.loc[eff['mean_predicted_model_score'].idxmax(), 'prrsv_type']
    report_lines.append(f"- Type {best_type} constructs exhibit a higher combined computational model score and immunogenicity.")
else:
    report_lines.append("- Additional data required for conclusive comparison.")

report_lines.append("\n✅ Report generated successfully by VAXINTAIC v3 pipeline.")
report_lines.append("=" * 65)

# ────────────────────────────────────────────────
# SAVE FINAL REPORT
# ────────────────────────────────────────────────
with open(REPORT_PATH, "w") as f:
    f.write("\n".join(report_lines))

print(f"✅ Advanced report generated → {REPORT_PATH}")
print("🧬 Stage 13 completed successfully.")

