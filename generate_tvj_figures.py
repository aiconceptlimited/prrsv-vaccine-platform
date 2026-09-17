#!/usr/bin/env python3

from pathlib import Path
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# VAXINTAIC — TVJ PUBLICATION FIGURE GENERATOR
# ============================================================

PROJECT_ROOT = Path.home() / "vaxintaic"
OUTPUT_DIR = PROJECT_ROOT / "VAXINTAIC_TVJ_Figures"

SCREENING_FILE = PROJECT_ROOT / "data/epitopes/epitope_predictions_v3.csv"
RANKING_FILE = PROJECT_ROOT / "data/candidates/top_epitopes_v3.csv"
CONSTRUCT_FILE = PROJECT_ROOT / "data/mrna/mrna_construct_metadata.csv"
CONSTRUCT_FASTA = PROJECT_ROOT / "data/mrna/mrna_constructs_v3.fasta"
NANOPARTICLE_FILE = PROJECT_ROOT / "data/models/nanoparticle_model_v3.csv"
MODEL_FILE = PROJECT_ROOT / "data/reports/model_scores_v3.csv"

ALIGN_I = PROJECT_ROOT / "data/alignments/prrsv_typei_aligned.fasta"
ALIGN_II = PROJECT_ROOT / "data/alignments/prrsv_typeii_aligned.fasta"

RUN_ID = "20260907T105902Z_7f5df939430f"


# ============================================================
# EXPECTED AUTHORITATIVE VALUES
# These are used ONLY for validation, not for plotting.
# ============================================================

EXPECTED_SCREENING = {
    "total": 14764,
    "Type I": 7332,
    "Type II": 7432,
}

EXPECTED_RANKING = {
    "Type I": 15,
    "Type II": 15,
}

EXPECTED_CONSTRUCT_LENGTH = 681


# ============================================================
# UTILITIES
# ============================================================

def fail(message):
    raise RuntimeError("\nERROR: " + message)


def require_file(path):
    if not path.exists():
        fail(f"Required authoritative file not found:\n{path}")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_figure(fig, stem):
    pdf = OUTPUT_DIR / f"{stem}.pdf"
    tiff = OUTPUT_DIR / f"{stem}.tiff"

    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(
        tiff,
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    return pdf, tiff


def read_fasta(path):
    records = []
    name = None
    seq = []

    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))

                name = line[1:].strip()
                seq = []
            else:
                seq.append(line)

    if name is not None:
        records.append((name, "".join(seq)))

    return records


def pairwise_alignment_statistics(records):
    if not records:
        fail("Alignment contains no sequences.")

    sequences = [seq.upper() for _, seq in records]

    lengths = {len(x) for x in sequences}

    if len(lengths) != 1:
        fail("Alignment sequences do not have identical lengths.")

    n = len(sequences)
    length = len(sequences[0])

    if n < 2:
        fail("At least two aligned sequences are required.")

    identities = []
    gaps = []

    for i in range(n):
        for j in range(i + 1, n):
            a = sequences[i]
            b = sequences[j]

            comparable = 0
            identical = 0

            for x, y in zip(a, b):
                if x == "-" or y == "-":
                    continue

                comparable += 1

                if x == y:
                    identical += 1

            if comparable:
                identities.append(identical / comparable)

            gap_fraction = (
                sum(x == "-" or y == "-" for x, y in zip(a, b))
                / length
            )

            gaps.append(gap_fraction)

    alignment_matrix = np.array(
        [
            [
                sum(
                    seqs_i[k] == seqs_j[k]
                    and seqs_i[k] != "-"
                    and seqs_j[k] != "-"
                    for k in range(length)
                )
                for seqs_j in sequences
            ]
            for seqs_i in sequences
        ],
        dtype=float,
    )

    conserved = 0
    variable = 0

    for pos in range(length):
        chars = [
            seq[pos]
            for seq in sequences
            if seq[pos] != "-"
        ]

        if not chars:
            continue

        if len(set(chars)) == 1:
            conserved += 1
        else:
            variable += 1

    # Dominance complement:
    # 1 - frequency of the most common residue at each position,
    # averaged across alignment positions.
    dominance_complements = []

    for pos in range(length):
        chars = [
            seq[pos]
            for seq in sequences
            if seq[pos] != "-"
        ]

        if not chars:
            continue

        counts = pd.Series(chars).value_counts(normalize=True)

        dominant_frequency = float(counts.iloc[0])
        dominance_complements.append(1.0 - dominant_frequency)

    diversity = (
        float(np.mean(dominance_complements))
        if dominance_complements
        else float("nan")
    )

    return {
        "n_sequences": n,
        "alignment_length": length,
        "mean_identity": float(np.mean(identities)),
        "mean_gap_fraction": float(np.mean(gaps)),
        "conserved_positions": conserved,
        "variable_positions": variable,
        "dominance_complement": diversity,
    }


def translate_dna(dna):
    table = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
        "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
        "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
        "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
        "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
        "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
        "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
        "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
        "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
        "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
        "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
        "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
        "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
        "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
        "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
    }

    if len(dna) % 3 != 0:
        fail("Construct nucleotide sequence length is not divisible by 3.")

    return "".join(
        table.get(dna[i:i+3], "X")
        for i in range(0, len(dna), 3)
    )


# ============================================================
# VALIDATION
# ============================================================

def validate_inputs():

    print("=" * 72)
    print("VAXINTAIC — TVJ PUBLICATION FIGURE GENERATOR")
    print("=" * 72)

    print(f"Project: {PROJECT_ROOT}")
    print(f"Output:  {OUTPUT_DIR}")
    print(f"Run ID:  {RUN_ID}")
    print()

    files = [
        SCREENING_FILE,
        RANKING_FILE,
        CONSTRUCT_FILE,
        CONSTRUCT_FASTA,
        NANOPARTICLE_FILE,
        MODEL_FILE,
        ALIGN_I,
        ALIGN_II,
    ]

    for path in files:
        require_file(path)
        print(f"✓ {path.relative_to(PROJECT_ROOT)}")

    print()

    screening = pd.read_csv(SCREENING_FILE)
    ranking = pd.read_csv(RANKING_FILE)
    construct = pd.read_csv(CONSTRUCT_FILE)
    nanoparticle = pd.read_csv(NANOPARTICLE_FILE)
    model = pd.read_csv(MODEL_FILE)

    required_screening = {
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
    }

    required_ranking = required_screening | {"VPI"}

    required_construct = {
        "construct_id",
        "prrsv_type",
        "num_epitopes",
        "length_nt",
        "GC_percent",
        "CAI",
    }

    required_nanoparticle = {
        "construct_id",
        "prrsv_type",
        "encapsulation_efficiency",
        "delivery_stability",
        "delivery_index",
    }

    required_model = {
        "construct_id",
        "prrsv_type",
        "mean_immunogenicity",
        "predicted_model_score",
        "vaccine_index",
    }

    for required, df, name in [
        (required_screening, screening, "screening"),
        (required_ranking, ranking, "ranking"),
        (required_construct, construct, "construct"),
        (required_nanoparticle, nanoparticle, "nanoparticle"),
        (required_model, model, "model"),
    ]:
        missing = required - set(df.columns)

        if missing:
            fail(
                f"{name} is missing required columns: "
                f"{sorted(missing)}"
            )

    # Screening count checks
    if len(screening) != EXPECTED_SCREENING["total"]:
        fail(
            f"Screening row count is {len(screening)}, "
            f"expected {EXPECTED_SCREENING['total']}."
        )

    counts = screening["prrsv_type"].value_counts().to_dict()

    for genotype, expected in EXPECTED_SCREENING.items():
        if genotype == "total":
            continue

        actual = int(counts.get(genotype, 0))

        if actual != expected:
            fail(
                f"{genotype} screening count is {actual}; "
                f"expected {expected}."
            )

    if not screening["sequence"].astype(str).str.len().eq(15).all():
        fail("Not all screening peptides are exactly 15 amino acids.")

    if not screening["sequence_type"].astype(str).eq("amino_acid").all():
        fail("Screening contains non-amino-acid records.")

    print(
        f"✓ Screening validated: "
        f"{len(screening):,} records "
        f"(Type I={counts['Type I']:,}; "
        f"Type II={counts['Type II']:,})"
    )

    # Ranking checks
    if len(ranking) != 30:
        fail(
            f"Ranking contains {len(ranking)} rows; "
            f"expected 30."
        )

    rank_counts = ranking["prrsv_type"].value_counts().to_dict()

    for genotype, expected in EXPECTED_RANKING.items():
        actual = int(rank_counts.get(genotype, 0))

        if actual != expected:
            fail(
                f"Ranking contains {actual} {genotype} candidates; "
                f"expected {expected}."
            )

    # Ensure ranked candidates exist in screening
    screening_keys = set(
        zip(
            screening["prrsv_type"],
            screening["start"],
            screening["end"],
            screening["sequence"],
        )
    )

    ranking_keys = set(
        zip(
            ranking["prrsv_type"],
            ranking["start"],
            ranking["end"],
            ranking["sequence"],
        )
    )

    missing = ranking_keys - screening_keys

    if missing:
        fail(
            f"{len(missing)} ranked candidates cannot be linked "
            f"to the authoritative screening output."
        )

    print("✓ Ranking validated: 30 records (15/type)")

    # Construct metadata
    if len(construct) != 2:
        fail(
            f"Construct metadata contains {len(construct)} rows; "
            f"expected 2."
        )

    if not construct["length_nt"].astype(int).eq(
        EXPECTED_CONSTRUCT_LENGTH
    ).all():
        fail("Construct metadata contains an unexpected nucleotide length.")

    print("✓ Construct metadata validated")

    # Cross-file identity
    construct_ids = set(construct["construct_id"])

    if set(model["construct_id"]) != construct_ids:
        fail("Model and construct metadata IDs do not match.")

    if set(nanoparticle["construct_id"]) != construct_ids:
        fail("Nanoparticle and construct metadata IDs do not match.")

    print("✓ Construct/model/nanoparticle linkage validated")

    return screening, ranking, construct, nanoparticle, model


# ============================================================
# CONSTRUCT FASTA VALIDATION
# ============================================================

def validate_construct_fasta(construct):

    records = read_fasta(CONSTRUCT_FASTA)

    if len(records) != 2:
        fail(
            f"Construct FASTA contains {len(records)} records; "
            f"expected 2."
        )

    fasta_lengths = {}

    for name, dna in records:

        dna = dna.upper().replace(" ", "")

        if len(dna) != EXPECTED_CONSTRUCT_LENGTH:
            fail(
                f"FASTA construct '{name}' has length {len(dna)} nt; "
                f"expected {EXPECTED_CONSTRUCT_LENGTH} nt."
            )

        if not dna.startswith("ATG"):
            fail(f"Construct '{name}' does not start with ATG.")

        if dna[-3:] not in {"TAA", "TAG", "TGA"}:
            fail(f"Construct '{name}' lacks a terminal stop codon.")

        protein = translate_dna(dna)

        if "*" in protein[:-1]:
            fail(
                f"Construct '{name}' contains an internal stop codon."
            )

        fasta_lengths[name] = len(dna)

    print("✓ Construct FASTA validated: 2 sequences, 681 nt each")
    print("✓ Terminal stop codons present")
    print("✓ Internal stop codons absent")

    return records


# ============================================================
# FIGURE 1
# ============================================================

def figure_1_workflow():

    stages = [
        "01  Fetch PRRSV sequences",
        "02  Multiple sequence alignment",
        "03  Sequence-derived peptide screening",
        "03b Epitope landscape visualization",
        "04  Candidate ranking and selection",
        "05  mRNA construct design",
        "06  In-silico nanoparticle modelling",
        "10  Immunogenicity scoring",
        "07  Construct-level computational prioritization",
        "09  Diversity analysis",
        "11  Sus scrofa codon adaptation",
        "13  Integrated report generation",
    ]

    fig, ax = plt.subplots(figsize=(8.2, 11))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(stages) + 1)
    ax.axis("off")

    ax.text(
        0.5,
        len(stages) + 0.55,
        "VAXINTAIC computational workflow",
        ha="center",
        va="center",
        fontsize=16,
        fontweight="bold",
    )

    for i, stage in enumerate(stages):
        y = len(stages) - i

        ax.text(
            0.5,
            y,
            stage,
            ha="center",
            va="center",
            fontsize=10,
            bbox=dict(
                boxstyle="round,pad=0.5",
                facecolor="white",
                edgecolor="black",
                linewidth=1.0,
            ),
        )

        if i < len(stages) - 1:
            ax.annotate(
                "",
                xy=(0.5, y - 0.72),
                xytext=(0.5, y - 0.28),
                arrowprops=dict(
                    arrowstyle="->",
                    linewidth=1.0,
                ),
            )

    ax.text(
        0.5,
        0.25,
        "Outputs are computational candidates and prioritization estimates "
        "requiring experimental assessment.",
        ha="center",
        va="center",
        fontsize=8.5,
    )

    return save_figure(fig, "Figure_1_VAXINTAIC_workflow")


# ============================================================
# FIGURE 2
# ============================================================

def figure_2_diversity(screening):

    stats_i = pairwise_alignment_statistics(read_fasta(ALIGN_I))
    stats_ii = pairwise_alignment_statistics(read_fasta(ALIGN_II))

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(9.2, 7.2),
    )

    genotypes = ["Type I", "Type II"]
    stats = [stats_i, stats_ii]

    identity_values = [
        stats_i["mean_identity"],
        stats_ii["mean_identity"],
    ]

    gap_values = [
        stats_i["mean_gap_fraction"],
        stats_ii["mean_gap_fraction"],
    ]

    conserved_values = [
        stats_i["conserved_positions"],
        stats_ii["conserved_positions"],
    ]

    variable_values = [
        stats_i["variable_positions"],
        stats_ii["variable_positions"],
    ]

    x = np.arange(2)

    axes[0, 0].bar(x, identity_values)
    axes[0, 0].set_xticks(x, genotypes)
    axes[0, 0].set_ylim(0, 1)
    axes[0, 0].set_ylabel("Mean pairwise identity")
    axes[0, 0].set_title("Sequence identity")

    axes[0, 1].bar(x, gap_values)
    axes[0, 1].set_xticks(x, genotypes)
    axes[0, 1].set_ylabel("Mean gap fraction")
    axes[0, 1].set_title("Alignment gaps")

    axes[1, 0].bar(
        x,
        conserved_values,
        label="Conserved",
    )
    axes[1, 0].bar(
        x,
        variable_values,
        bottom=conserved_values,
        label="Variable",
    )
    axes[1, 0].set_xticks(x, genotypes)
    axes[1, 0].set_ylabel("Alignment positions")
    axes[1, 0].set_title("Conserved and variable positions")
    axes[1, 0].legend(frameon=False)

    diversity_values = [
        stats_i["dominance_complement"],
        stats_ii["dominance_complement"],
    ]

    axes[1, 1].bar(x, diversity_values)
    axes[1, 1].set_xticks(x, genotypes)
    axes[1, 1].set_ylabel("Dominance-complement")
    axes[1, 1].set_title("Sequence diversity")

    fig.suptitle(
        "PRRSV sequence diversity and alignment characteristics",
        fontsize=14,
        fontweight="bold",
    )

    fig.tight_layout()

    return save_figure(fig, "Figure_2_PRSV_diversity_alignment")


# ============================================================
# FIGURE 3
# ============================================================

def figure_3_candidate_landscape(screening, ranking):

    fig, ax = plt.subplots(figsize=(10, 5.8))

    for genotype in ["Type I", "Type II"]:

        subset = screening[
            screening["prrsv_type"] == genotype
        ].copy()

        ax.scatter(
            subset["start"],
            subset["epitope_score"],
            s=10,
            alpha=0.30,
            label=f"{genotype} screening windows",
        )

        ranked = ranking[
            ranking["prrsv_type"] == genotype
        ].copy()

        # Match by stable biological key rather than dataframe index.
        ranked_keys = set(
            zip(
                ranked["start"],
                ranked["end"],
                ranked["sequence"],
            )
        )

        selected = subset[
            subset.apply(
                lambda row: (
                    row["start"],
                    row["end"],
                    row["sequence"],
                ) in ranked_keys,
                axis=1,
            )
        ]

        ax.scatter(
            selected["start"],
            selected["epitope_score"],
            s=42,
            edgecolors="black",
            linewidths=0.7,
            label=f"{genotype} selected candidates",
        )

    ax.set_xlabel("Peptide start position")
    ax.set_ylabel("Epitope score")
    ax.set_title(
        "Sequence-derived peptide candidate landscape",
        fontweight="bold",
    )

    ax.legend(
        frameon=False,
        fontsize=8,
    )

    fig.tight_layout()

    return save_figure(fig, "Figure_3_peptide_candidate_landscape")


# ============================================================
# FIGURE 4
# ============================================================

def figure_4_prioritization(construct, nanoparticle, model):

    merged = construct.merge(
        nanoparticle[
            [
                "construct_id",
                "encapsulation_efficiency",
                "delivery_stability",
                "delivery_index",
            ]
        ],
        on="construct_id",
        validate="one_to_one",
    )

    merged = merged.merge(
        model[
            [
                "construct_id",
                "mean_immunogenicity",
                "predicted_model_score",
                "vaccine_index",
            ]
        ],
        on="construct_id",
        validate="one_to_one",
    )

    merged = merged.sort_values("prrsv_type")

    labels = merged["prrsv_type"].tolist()

    x = np.arange(len(merged))
    width = 0.18

    fig, ax = plt.subplots(figsize=(10, 6))

    metrics = [
        (
            "encapsulation_efficiency",
            "Encapsulation efficiency",
        ),
        (
            "delivery_stability",
            "Delivery stability",
        ),
        (
            "mean_immunogenicity",
            "Mean immunogenicity",
        ),
        (
            "predicted_model_score",
            "Computational model score",
        ),
        (
            "vaccine_index",
            "Integrated prioritization index",
        ),
    ]

    for i, (column, label) in enumerate(metrics):
        ax.bar(
            x + (i - 2) * width,
            merged[column].astype(float),
            width,
            label=label,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Modelled / normalized value")
    ax.set_title(
        "Construct-level computational prioritization",
        fontweight="bold",
    )

    ax.legend(
        frameon=False,
        fontsize=8,
        ncol=2,
    )

    fig.tight_layout()

    return save_figure(fig, "Figure_4_construct_prioritization")


# ============================================================
# MANIFEST
# ============================================================

def create_manifest(files, figures):

    manifest = {
        "project": "VAXINTAIC",
        "purpose": "The Veterinary Journal publication figure generation",
        "pipeline_run_id": RUN_ID,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(PROJECT_ROOT),
        "output_directory": str(OUTPUT_DIR),
        "python_version": sys.version,
        "platform": platform.platform(),
        "random_scientific_data_generated": False,
        "inputs": {},
        "outputs": {},
    }

    for path in files:
        manifest["inputs"][str(path.relative_to(PROJECT_ROOT))] = {
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }

    for path in figures:
        manifest["outputs"][path.name] = {
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }

    manifest_path = OUTPUT_DIR / "MANIFEST.json"

    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    return manifest_path


# ============================================================
# VALIDATION REPORT
# ============================================================

def write_validation_report(screening, ranking, construct):

    stats_i = pairwise_alignment_statistics(read_fasta(ALIGN_I))
    stats_ii = pairwise_alignment_statistics(read_fasta(ALIGN_II))

    lines = [
        "VAXINTAIC TVJ PUBLICATION FIGURE VALIDATION",
        "=" * 55,
        "",
        f"Pipeline run ID: {RUN_ID}",
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}",
        "",
        "SOURCE DATA CHECKS",
        "------------------",
        f"Screening records: {len(screening):,} ✓",
        f"Type I screening records: "
        f"{(screening.prrsv_type == 'Type I').sum():,} ✓",
        f"Type II screening records: "
        f"{(screening.prrsv_type == 'Type II').sum():,} ✓",
        f"Ranked candidates: {len(ranking):,} ✓",
        f"Type I ranked candidates: "
        f"{(ranking.prrsv_type == 'Type I').sum():,} ✓",
        f"Type II ranked candidates: "
        f"{(ranking.prrsv_type == 'Type II').sum():,} ✓",
        "",
        "SEQUENCE CHECKS",
        "---------------",
        "All screening peptides are 15 aa: PASS ✓",
        "All screening records are amino-acid records: PASS ✓",
        "Ranked candidates linked to screening data: PASS ✓",
        "Construct metadata records: 2 ✓",
        "Construct length: 681 nt each ✓",
        "Construct FASTA terminal stop codons: PASS ✓",
        "Construct FASTA internal stop codons: 0 ✓",
        "",
        "ALIGNMENT-DERIVED VALUES",
        "-------------------------",
        f"Type I mean identity: "
        f"{stats_i['mean_identity']:.6f}",
        f"Type I mean gap fraction: "
        f"{stats_i['mean_gap_fraction']:.6f}",
        f"Type I conserved positions: "
        f"{stats_i['conserved_positions']}",
        f"Type I variable positions: "
        f"{stats_i['variable_positions']}",
        f"Type I dominance-complement: "
        f"{stats_i['dominance_complement']:.6f}",
        "",
        f"Type II mean identity: "
        f"{stats_ii['mean_identity']:.6f}",
        f"Type II mean gap fraction: "
        f"{stats_ii['mean_gap_fraction']:.6f}",
        f"Type II conserved positions: "
        f"{stats_ii['conserved_positions']}",
        f"Type II variable positions: "
        f"{stats_ii['variable_positions']}",
        f"Type II dominance-complement: "
        f"{stats_ii['dominance_complement']:.6f}",
        "",
        "FIGURE GENERATION",
        "------------------",
        "Figure 1: PASS ✓",
        "Figure 2: PASS ✓",
        "Figure 3: PASS ✓",
        "Figure 4: PASS ✓",
        "",
        "Synthetic/random scientific data generated: NO ✓",
        "",
        "STATUS: AUTHORITATIVE SOURCE FILES VALIDATED",
    ]

    report = OUTPUT_DIR / "VALIDATION_REPORT.txt"

    with open(report, "w") as fh:
        fh.write("\n".join(lines))

    return report


# ============================================================
# MAIN
# ============================================================

def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    (
        screening,
        ranking,
        construct,
        nanoparticle,
        model,
    ) = validate_inputs()

    validate_construct_fasta(construct)

    print()
    print("Generating figures...")

    figures = []

    for result in [
        figure_1_workflow(),
        figure_2_diversity(screening),
        figure_3_candidate_landscape(
            screening,
            ranking,
        ),
        figure_4_prioritization(
            construct,
            nanoparticle,
            model,
        ),
    ]:
        figures.extend(result)

    input_files = [
        SCREENING_FILE,
        RANKING_FILE,
        CONSTRUCT_FILE,
        CONSTRUCT_FASTA,
        NANOPARTICLE_FILE,
        MODEL_FILE,
        ALIGN_I,
        ALIGN_II,
    ]

    report = write_validation_report(
        screening,
        ranking,
        construct,
    )

    manifest = create_manifest(
        input_files,
        figures,
    )

    print()
    print("=" * 72)
    print("FIGURE GENERATION COMPLETE")
    print("=" * 72)

    for path in figures:
        print(f"✓ {path}")

    print(f"✓ {report}")
    print(f"✓ {manifest}")

    print()
    print("No synthetic/random scientific data were generated.")
    print("STATUS: AUTHORITATIVE SOURCE FILES VALIDATED")


if __name__ == "__main__":
    main()

