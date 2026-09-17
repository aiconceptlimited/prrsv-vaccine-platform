#!/usr/bin/env python3

from pathlib import Path
import hashlib
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = Path.home() / "vaxintaic"
OUT_DIR = BASE_DIR / "VAXINTAIC_TVJ_Figures"

ALIGNMENT_SUMMARY = BASE_DIR / "data/alignments/alignment_summary_v3.csv"
DIVERSITY_SUMMARY = BASE_DIR / "data/analysis/diversity_summary_v3.csv"
TYPE_I_PROFILE = BASE_DIR / "data/analysis/diversity_profile_typei_v3.csv"
TYPE_II_PROFILE = BASE_DIR / "data/analysis/diversity_profile_typeii_v3.csv"

PDF_OUT = OUT_DIR / "Figure_2_PRSV_diversity_alignment.pdf"
TIFF_OUT = OUT_DIR / "Figure_2_PRSV_diversity_alignment.tiff"
VALIDATION_OUT = OUT_DIR / "FIGURE_2_VALIDATION.txt"

EXPECTED_TYPES = ["Type I", "Type II"]

EXPECTED_ALIGNMENT = {
    "Type I": {
        "num_sequences": 40,
        "alignment_length": 619,
        "gap_fraction": 0.0250,
        "avg_identity": 0.7791,
    },
    "Type II": {
        "num_sequences": 40,
        "alignment_length": 603,
        "gap_fraction": 0.0001,
        "avg_identity": 0.9335,
    },
}

EXPECTED_DIVERSITY = {
    "Type I": {
        "mean_diversity": 0.1618,
        "conserved_positions": 321,
        "variable_positions": 146,
    },
    "Type II": {
        "mean_diversity": 0.0395,
        "conserved_positions": 520,
        "variable_positions": 17,
    },
}

TOL = 1e-9


def fail(message):
    raise RuntimeError("FAIL-CLOSED: " + message)


def require_file(path):
    if not path.exists():
        fail(f"Missing authoritative file: {path}")
    if path.stat().st_size == 0:
        fail(f"Empty authoritative file: {path}")


def require_columns(df, columns, label):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        fail(f"{label} missing columns: {missing}")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def check_close(actual, expected, label):
    if not np.isclose(float(actual), float(expected), rtol=0, atol=TOL):
        fail(f"{label}: {actual} != {expected}")


# ------------------------------------------------------------
# Stage 02 — authoritative alignment summary
# ------------------------------------------------------------

require_file(ALIGNMENT_SUMMARY)

alignment = pd.read_csv(ALIGNMENT_SUMMARY)

require_columns(
    alignment,
    [
        "type",
        "num_sequences",
        "alignment_length",
        "gap_fraction",
        "avg_identity",
    ],
    "Stage 02 alignment summary",
)

if alignment["type"].duplicated().any():
    fail("Duplicate genotype rows in Stage 02 summary")

if set(alignment["type"]) != set(EXPECTED_TYPES):
    fail(f"Unexpected genotypes in Stage 02 summary: {alignment['type'].tolist()}")

alignment = alignment.set_index("type").loc[EXPECTED_TYPES]

for t in EXPECTED_TYPES:
    row = alignment.loc[t]
    exp = EXPECTED_ALIGNMENT[t]

    if int(row["num_sequences"]) != exp["num_sequences"]:
        fail(f"{t} sequence count mismatch")

    if int(row["alignment_length"]) != exp["alignment_length"]:
        fail(f"{t} alignment length mismatch")

    check_close(row["gap_fraction"], exp["gap_fraction"], f"{t} gap fraction")
    check_close(row["avg_identity"], exp["avg_identity"], f"{t} mean identity")


# ------------------------------------------------------------
# Stage 09 — authoritative diversity summary
# ------------------------------------------------------------

require_file(DIVERSITY_SUMMARY)

diversity = pd.read_csv(DIVERSITY_SUMMARY)

require_columns(
    diversity,
    [
        "prrsv_type",
        "mean_diversity",
        "conserved_positions",
        "variable_positions",
    ],
    "Stage 09 diversity summary",
)

if diversity["prrsv_type"].duplicated().any():
    fail("Duplicate genotype rows in Stage 09 summary")

if set(diversity["prrsv_type"]) != set(EXPECTED_TYPES):
    fail(
        f"Unexpected genotypes in Stage 09 summary: "
        f"{diversity['prrsv_type'].tolist()}"
    )

diversity = diversity.set_index("prrsv_type").loc[EXPECTED_TYPES]

for t in EXPECTED_TYPES:
    row = diversity.loc[t]
    exp = EXPECTED_DIVERSITY[t]

    check_close(row["mean_diversity"], exp["mean_diversity"], f"{t} mean diversity")

    if int(row["conserved_positions"]) != exp["conserved_positions"]:
        fail(f"{t} conserved-position count mismatch")

    if int(row["variable_positions"]) != exp["variable_positions"]:
        fail(f"{t} variable-position count mismatch")


# ------------------------------------------------------------
# Stage 09 — authoritative diversity profiles
# ------------------------------------------------------------

require_file(TYPE_I_PROFILE)
require_file(TYPE_II_PROFILE)

profiles = {
    "Type I": pd.read_csv(TYPE_I_PROFILE),
    "Type II": pd.read_csv(TYPE_II_PROFILE),
}

for t, df in profiles.items():

    require_columns(
        df,
        ["position", "diversity"],
        f"{t} diversity profile",
    )

    if df.empty:
        fail(f"{t} diversity profile is empty")

    if df["position"].duplicated().any():
        fail(f"{t} diversity profile has duplicate positions")

    if df["position"].isna().any() or df["diversity"].isna().any():
        fail(f"{t} diversity profile contains missing values")

    values = df["diversity"].to_numpy(dtype=float)

    if not np.isfinite(values).all():
        fail(f"{t} diversity profile contains non-finite values")

    if (values < 0).any() or (values > 1).any():
        fail(f"{t} diversity values outside 0–1")

    positions = df["position"].astype(int).to_numpy()

    if not np.array_equal(positions, np.arange(1, len(df) + 1)):
        fail(f"{t} profile positions are not continuous 1-based positions")

    # Verify that Stage 09 summary and profile agree.
    mean_profile = df["diversity"].mean()
    conserved_profile = int((df["diversity"] < 0.1).sum())
    variable_profile = int((df["diversity"] > 0.3).sum())

    # Stage 09 writes mean_diversity rounded to 4 decimal places.
    # Therefore compare the full profile mean to the published
    # Stage 09 summary after rounding to the same precision.
    summary_mean = round(float(diversity.loc[t, "mean_diversity"]), 4)
    profile_mean_rounded = round(float(mean_profile), 4)

    if profile_mean_rounded != summary_mean:
        fail(
            f"{t} profile mean vs summary mismatch after Stage 09 "
            f"rounding: {profile_mean_rounded} != {summary_mean}"
        )

    if conserved_profile != int(diversity.loc[t, "conserved_positions"]):
        fail(f"{t} conserved count differs between profile and summary")

    if variable_profile != int(diversity.loc[t, "variable_positions"]):
        fail(f"{t} variable count differs between profile and summary")


# ------------------------------------------------------------
# Generate Figure 2
# ------------------------------------------------------------

OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

fig = plt.figure(figsize=(11, 8.2))

gs = fig.add_gridspec(
    2,
    2,
    height_ratios=[1.0, 1.4],
    hspace=0.42,
    wspace=0.32,
)

# ------------------------------------------------------------
# A — Alignment characteristics
# ------------------------------------------------------------

ax1 = fig.add_subplot(gs[0, 0])

x = np.arange(2)
width = 0.36

identity = [
    float(alignment.loc["Type I", "avg_identity"]),
    float(alignment.loc["Type II", "avg_identity"]),
]

gap = [
    float(alignment.loc["Type I", "gap_fraction"]),
    float(alignment.loc["Type II", "gap_fraction"]),
]

ax1.bar(
    x - width / 2,
    identity,
    width,
    label="Mean pairwise identity",
)

ax1.bar(
    x + width / 2,
    gap,
    width,
    label="Gap fraction",
)

ax1.set_xticks(x)
ax1.set_xticklabels(["Type I", "Type II"])
ax1.set_ylabel("Fraction")
ax1.set_ylim(0, 1.08)
ax1.set_title("A. Alignment characteristics")
ax1.legend(frameon=False, fontsize=8)

for i, v in enumerate(identity):
    ax1.text(
        i - width / 2,
        v + 0.025,
        f"{v:.4f}",
        ha="center",
        va="bottom",
        fontsize=8,
    )

for i, v in enumerate(gap):
    ax1.text(
        i + width / 2,
        v + 0.025,
        f"{v:.4f}",
        ha="center",
        va="bottom",
        fontsize=8,
    )


# ------------------------------------------------------------
# B — Diversity summary
# ------------------------------------------------------------

ax2 = fig.add_subplot(gs[0, 1])

mean_div = [
    float(diversity.loc["Type I", "mean_diversity"]),
    float(diversity.loc["Type II", "mean_diversity"]),
]

ax2.bar(x, mean_div, width=0.55)

ax2.set_xticks(x)
ax2.set_xticklabels(["Type I", "Type II"])
ax2.set_ylabel("Mean diversity")
ax2.set_title("B. Dominance-complement diversity")

ymax = max(mean_div) * 1.45
ax2.set_ylim(0, ymax)

for i, v in enumerate(mean_div):
    ax2.text(
        i,
        v + ymax * 0.035,
        f"{v:.4f}",
        ha="center",
        va="bottom",
        fontsize=9,
    )

for i, t in enumerate(EXPECTED_TYPES):
    conserved = int(diversity.loc[t, "conserved_positions"])
    variable = int(diversity.loc[t, "variable_positions"])

    ax2.text(
        i,
        ymax * 0.93,
        f"Conserved: {conserved}\nVariable: {variable}",
        ha="center",
        va="top",
        fontsize=8,
    )


# ------------------------------------------------------------
# C — Stage 09 per-position profiles
# ------------------------------------------------------------

ax3 = fig.add_subplot(gs[1, :])

for t in EXPECTED_TYPES:
    df = profiles[t]

    ax3.plot(
        df["position"],
        df["diversity"],
        linewidth=1.15,
        label=t,
    )

ax3.axhline(
    0.1,
    linestyle="--",
    linewidth=0.8,
    label="Conserved threshold (0.1)",
)

ax3.axhline(
    0.3,
    linestyle=":",
    linewidth=0.8,
    label="Variable threshold (0.3)",
)

ax3.set_xlabel("Alignment position")
ax3.set_ylabel("Dominance-complement diversity")
ax3.set_title("C. Per-position diversity profiles")

profile_max = max(
    profiles["Type I"]["diversity"].max(),
    profiles["Type II"]["diversity"].max(),
    0.3,
)

ax3.set_ylim(0, profile_max * 1.08)

ax3.legend(
    frameon=False,
    ncol=2,
    fontsize=8,
    loc="upper right",
)

ax3.text(
    0.01,
    -0.18,
    "Metric = 1 − dominant_frequency. "
    "Conserved: diversity < 0.1; "
    "variable: diversity > 0.3.",
    transform=ax3.transAxes,
    fontsize=8,
    va="top",
)

fig.suptitle(
    "PRRSV sequence diversity and alignment characteristics",
    fontsize=13,
    y=0.98,
)

fig.text(
    0.5,
    0.015,
    "VAXINTAIC — authoritative Stage 02 and Stage 09 outputs",
    ha="center",
    fontsize=8,
)

fig.savefig(
    PDF_OUT,
    format="pdf",
    bbox_inches="tight",
)

fig.savefig(
    TIFF_OUT,
    format="tiff",
    dpi=600,
    bbox_inches="tight",
)

plt.close(fig)


# ------------------------------------------------------------
# Validation report
# ------------------------------------------------------------

now = datetime.datetime.now(datetime.timezone.utc).isoformat()

lines = [
    "VAXINTAIC TVJ FIGURE 2 VALIDATION",
    "=" * 60,
    "",
    f"Generated UTC: {now}",
    "",
    "AUTHORITATIVE INPUT FILES",
    "-------------------------",
]

for p in [
    ALIGNMENT_SUMMARY,
    DIVERSITY_SUMMARY,
    TYPE_I_PROFILE,
    TYPE_II_PROFILE,
]:
    lines.append(f"{p}: PASS")
    lines.append(f"SHA-256: {sha256_file(p)}")

lines += [
    "",
    "STAGE 02 ALIGNMENT VALUES",
    "--------------------------",
]

for t in EXPECTED_TYPES:
    r = alignment.loc[t]
    lines.append(
        f"{t}: n={int(r['num_sequences'])}, "
        f"alignment_length={int(r['alignment_length'])}, "
        f"gap_fraction={float(r['gap_fraction']):.4f}, "
        f"mean_identity={float(r['avg_identity']):.4f}"
    )

lines += [
    "",
    "STAGE 09 DIVERSITY VALUES",
    "-------------------------",
]

for t in EXPECTED_TYPES:
    r = diversity.loc[t]
    lines.append(
        f"{t}: mean_diversity={float(r['mean_diversity']):.4f}, "
        f"conserved_positions={int(r['conserved_positions'])}, "
        f"variable_positions={int(r['variable_positions'])}"
    )

lines += [
    "",
    "DIVERSITY DEFINITION",
    "--------------------",
    "Metric: dominance_complement",
    "Formula: 1 - dominant_frequency",
    "Conserved threshold: diversity < 0.1",
    "Variable threshold: diversity > 0.3",
    "",
    "SYNTHETIC/RANDOM SCIENTIFIC DATA",
    "--------------------------------",
    "NO — all plotted scientific values originate from",
    "authoritative Stage 02/Stage 09 output files.",
    "",
    "OUTPUTS",
    "-------",
    f"PDF: {PDF_OUT}",
    f"TIFF: {TIFF_OUT}",
    "",
    "STATUS: PASS — FIGURE 2 IS SOURCE-FAITHFUL",
]

VALIDATION_OUT.write_text(
    "\n".join(lines) + "\n",
    encoding="utf-8",
)

print("")
print("==============================================")
print(" FIGURE 2 GENERATED SUCCESSFULLY")
print("==============================================")
print(f"✓ {PDF_OUT}")
print(f"✓ {TIFF_OUT}")
print(f"✓ {VALIDATION_OUT}")
print("")
print("AUTHORITATIVE VALUES:")
print("---------------------")

for t in EXPECTED_TYPES:
    a = alignment.loc[t]
    d = diversity.loc[t]

    print(
        f"{t}: identity={a['avg_identity']:.4f}, "
        f"gap={a['gap_fraction']:.4f}, "
        f"diversity={d['mean_diversity']:.4f}, "
        f"conserved={int(d['conserved_positions'])}, "
        f"variable={int(d['variable_positions'])}"
    )

print("")
print("STATUS: PASS — NO SYNTHETIC SCIENTIFIC DATA")
