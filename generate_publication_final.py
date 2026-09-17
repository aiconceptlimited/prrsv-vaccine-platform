#!/usr/bin/env python3
"""
VAXINTAIC TVJ Publication Generator (Production Version)
Run from: ~/vaxintaic
Output: ~/vaxintaic/VAXINTAIC_TVJ_Submission/
"""

import os
import sys
import json
import hashlib
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = PROJECT_ROOT / "VAXINTAIC_TVJ_Submission"

# AUTHORITATIVE FILE PATHS (actual VPS locations)
FILES = {
    'screening': PROJECT_ROOT / "data/epitopes/epitope_predictions_v3.csv",
    'ranking': PROJECT_ROOT / "data/candidates/top_epitopes_v3.csv",
    'constructs': PROJECT_ROOT / "data/mrna/mrna_construct_metadata.csv",
    'nanoparticle': PROJECT_ROOT / "data/models/nanoparticle_model_v3.csv",
    'immunogenicity': PROJECT_ROOT / "data/epitopes/immunogenicity_scored.csv",
    'model': PROJECT_ROOT / "data/reports/model_scores_v3.csv",
    'alignment': PROJECT_ROOT / "data/alignments/alignment_summary_v3.csv",
    'diversity': PROJECT_ROOT / "data/analysis/diversity_summary_v3.csv",
    'cai': PROJECT_ROOT / "data/mrna/codon_adaptation_v3.csv",
    'report': PROJECT_ROOT / "data/validation_report.json",
}

# Expected summary values (for cross-validation)
EXPECTED = {
    'screening_total': 14764,
    'screening_type_i': 7332,
    'screening_type_ii': 7432,
    'ranking_type_i': 15,
    'ranking_type_ii': 15,
    'construct_length': 681,
}

# ============================================================================
# VALIDATION
# ============================================================================

def validate_file(path, desc):
    if not path.exists():
        raise FileNotFoundError(f"ERROR: {desc} not found at {path}")
    print(f"✓ {desc}: {path}")

def sha256_dataframe(df):
    return hashlib.sha256(df.to_csv(index=False).encode()).hexdigest()[:16]

# ============================================================================
# LOAD DATA
# ============================================================================

def load_data():
    data = {}
    hashes = {}

    validate_file(FILES['screening'], "Screening")
    data['screening'] = pd.read_csv(FILES['screening'])
    hashes['screening'] = sha256_dataframe(data['screening'])
    print(f"  Screening: {len(data['screening'])} records")

    validate_file(FILES['ranking'], "Ranking")
    data['ranking'] = pd.read_csv(FILES['ranking'])
    hashes['ranking'] = sha256_dataframe(data['ranking'])
    print(f"  Ranking: {len(data['ranking'])} records")

    validate_file(FILES['constructs'], "Constructs")
    data['constructs'] = pd.read_csv(FILES['constructs'])
    hashes['constructs'] = sha256_dataframe(data['constructs'])
    print(f"  Constructs: {len(data['constructs'])} records")

    validate_file(FILES['nanoparticle'], "Nanoparticle")
    data['nanoparticle'] = pd.read_csv(FILES['nanoparticle'])
    hashes['nanoparticle'] = sha256_dataframe(data['nanoparticle'])

    validate_file(FILES['immunogenicity'], "Immunogenicity")
    data['immunogenicity'] = pd.read_csv(FILES['immunogenicity'])
    hashes['immunogenicity'] = sha256_dataframe(data['immunogenicity'])

    validate_file(FILES['model'], "Model")
    data['model'] = pd.read_csv(FILES['model'])
    hashes['model'] = sha256_dataframe(data['model'])

    validate_file(FILES['diversity'], "Diversity")
    data['diversity'] = pd.read_csv(FILES['diversity'])
    hashes['diversity'] = sha256_dataframe(data['diversity'])

    validate_file(FILES['cai'], "CAI")
    data['cai'] = pd.read_csv(FILES['cai'])
    hashes['cai'] = sha256_dataframe(data['cai'])

    if FILES['report'].exists():
        with open(FILES['report'], 'r') as f:
            data['report'] = json.load(f)
        hashes['report'] = hashlib.sha256(open(FILES['report'], 'rb').read()).hexdigest()[:16]
    else:
        data['report'] = {}
        hashes['report'] = None

    return data, hashes

# ============================================================================
# VERIFY DATA
# ============================================================================

def verify_data(data):
    print('\n' + '=' * 70)
    print('VERIFYING DATA QUALITY')
    print('=' * 70 + '\n')

    errors = []

    # 1. Screening counts
    type_i = data['screening'][data['screening']['prrsv_type'] == 'Type I']
    type_ii = data['screening'][data['screening']['prrsv_type'] == 'Type II']

    if len(data['screening']) != EXPECTED['screening_total']:
        errors.append(f"Screening total: expected {EXPECTED['screening_total']}, got {len(data['screening'])}")
    if len(type_i) != EXPECTED['screening_type_i']:
        errors.append(f"Type I screening: expected {EXPECTED['screening_type_i']}, got {len(type_i)}")
    if len(type_ii) != EXPECTED['screening_type_ii']:
        errors.append(f"Type II screening: expected {EXPECTED['screening_type_ii']}, got {len(type_ii)}")

    # 2. Peptide length
    invalid = data['screening'][data['screening']['sequence'].str.len() != 15]
    if len(invalid) > 0:
        errors.append(f"{len(invalid)} peptides not 15 aa")

    # 3. Ranking counts
    type_i_rank = data['ranking'][data['ranking']['prrsv_type'] == 'Type I']
    type_ii_rank = data['ranking'][data['ranking']['prrsv_type'] == 'Type II']

    if len(type_i_rank) != EXPECTED['ranking_type_i']:
        errors.append(f"Type I ranking: expected {EXPECTED['ranking_type_i']}, got {len(type_i_rank)}")
    if len(type_ii_rank) != EXPECTED['ranking_type_ii']:
        errors.append(f"Type II ranking: expected {EXPECTED['ranking_type_ii']}, got {len(type_ii_rank)}")

    # 4. Verify ranked candidates exist in screening
    screening_keys = set(data['screening']['sequence'] + '_' + data['screening']['start'].astype(str))
    rank_keys = set(data['ranking']['sequence'] + '_' + data['ranking']['start'].astype(str))
    missing = rank_keys - screening_keys
    if missing:
        errors.append(f"{len(missing)} ranked candidates not found in screening")

    # 5. Construct validation
    for _, row in data['constructs'].iterrows():
        nt_len = int(row.get('length_nt', 0))
        if nt_len != EXPECTED['construct_length']:
            errors.append(f"Construct {row.get('construct_id')}: expected {EXPECTED['construct_length']} nt, got {nt_len}")

    if errors:
        print('!' * 70)
        for err in errors:
            print(f'  ✗ {err}')
        raise ValueError(f'{len(errors)} validation errors')

    print('✓ All data validation checks passed')
    print(f'  Screening: {len(data["screening"])} records')
    print(f'  Type I: {len(type_i)}, Type II: {len(type_ii)}')
    print(f'  Top candidates: Type I {len(type_i_rank)}, Type II {len(type_ii_rank)}')

# ============================================================================
# GENERATE FIGURES
# ============================================================================

def set_style():
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = 11
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.bbox'] = 'tight'
    plt.rcParams['pdf.fonttype'] = 42

def figure_1_workflow():
    set_style()
    fig, ax = plt.subplots(1, 1, figsize=(14, 11))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 13)
    ax.axis('off')

    stages = [
        ('01', 'Sequence Acquisition', '80 PRRSV sequences\n(40 Type I, 40 Type II)', '#2E86AB'),
        ('02', 'Multiple Sequence Alignment', 'MAFFT v7', '#2E86AB'),
        ('03', 'Peptide Screening', 'Hydrophilicity, flexibility,\nconservation', '#2E86AB'),
        ('03b', 'Epitope Landscape Plotting', 'Visualization of\ncandidate scores', '#2E86AB'),
        ('04', 'Candidate Ranking', 'Top 15 candidates per\ngenotype', '#2E86AB'),
        ('05', 'mRNA Construct Design', '681 nt constructs\n(Type I, Type II)', '#F4A261'),
        ('06', 'Nanoparticle Delivery\nModelling', 'In-silico encapsulation\nand stability', '#E76F51'),
        ('10', 'Immunogenicity Scoring', 'VAXINTAIC computational\nscores', '#E76F51'),
        ('07', 'Construct-Level\nPrioritization', 'GradientBoostingRegressor', '#E76F51'),
        ('09', 'Diversity Analysis', 'Dominance-complement\nmetric', '#F4A261'),
        ('11', 'Codon Adaptation\nAnalysis', 'Sus scrofa CAI', '#F4A261'),
        ('13', 'Automated Reporting', 'Integrated candidate\ncomparison', '#2E9B5E'),
    ]

    x_pos = [2, 7, 12]
    y_pos = [11, 9, 7, 5, 3, 1]

    for i, (sid, name, desc, color) in enumerate(stages):
        col, row = i % 3, i // 3
        if row >= len(y_pos): continue
        x, y = x_pos[col], y_pos[row]
        box = FancyBboxPatch((x - 1.8, y - 0.7), 3.6, 1.4,
                             boxstyle='round,pad=0.1', facecolor=color,
                             edgecolor='black', linewidth=1.5, alpha=0.85)
        ax.add_patch(box)
        ax.text(x, y + 0.35, f'{sid}:', ha='center', va='center', fontweight='bold', fontsize=8)
        ax.text(x, y - 0.05, name, ha='center', va='center', fontweight='bold', fontsize=8)
        ax.text(x, y - 0.45, desc, ha='center', va='center', fontsize=7)

    # Clean orthogonal workflow connectors.
    # Horizontal transitions connect adjacent boxes directly.
    # Row transitions use a right-angle path rather than diagonal arrows.
    for i in range(len(stages) - 1):
        col, row = i % 3, i // 3
        nc, nr = (i + 1) % 3, (i + 1) // 3

        x1, y1 = x_pos[col], y_pos[row]
        x2, y2 = x_pos[nc], y_pos[nr]

        if nr == row:
            ax.add_patch(
                FancyArrowPatch(
                    (x1 + 1.8, y1),
                    (x2 - 1.8, y2),
                    arrowstyle='->',
                    mutation_scale=13,
                    color='gray',
                    linewidth=1.4
                )
            )
        else:
            y_mid = (y1 + y2) / 2

            ax.plot(
                [x1, x1],
                [y1 - 0.7, y_mid],
                color='gray',
                linewidth=1.4
            )

            ax.plot(
                [x1, x2],
                [y_mid, y_mid],
                color='gray',
                linewidth=1.4
            )

            ax.add_patch(
                FancyArrowPatch(
                    (x2, y_mid),
                    (x2, y2 + 0.7),
                    arrowstyle='->',
                    mutation_scale=13,
                    color='gray',
                    linewidth=1.4
                )
            )

    legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor='#2E86AB', alpha=0.85, label='Sequence Analysis'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#F4A261', alpha=0.85, label='Construct Design'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#E76F51', alpha=0.85, label='Computational Characterization'),
        plt.Rectangle((0, 0), 1, 1, facecolor='#2E9B5E', alpha=0.85, label='Output'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', frameon=True, fontsize=9)
    ax.text(7, 12.6, 'VAXINTAIC Computational Workflow', ha='center', va='center', fontsize=14, fontweight='bold')
    ax.text(7, 0.3, 'Outputs are computational prioritization measures requiring experimental validation',
            ha='center', va='center', fontsize=9, style='italic')

    # Fixed publication framing for stable PDF/TIFF rendering.
    fig = plt.gcf()
    fig.subplots_adjust(
        left=0.055,
        right=0.945,
        top=0.91,
        bottom=0.10
    )

    (OUTPUT_ROOT / 'figures').mkdir(parents=True, exist_ok=True)
    for ext in ['pdf', 'tiff']:
        plt.savefig(
            OUTPUT_ROOT / 'figures' / f'Figure1_workflow.{ext}',
            dpi=300,
            bbox_inches=None
        )
    print('✓ Figure 1 saved')
    plt.close()

def figure_2_diversity():
    set_style()

    # Publication layout: figure-level heading separated from panel labels.
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.2))

    # Authoritative VAXINTAIC analysis outputs
    alignment = pd.read_csv("data/alignments/alignment_summary_v3.csv")
    diversity_summary = pd.read_csv("data/analysis/diversity_summary_v3.csv")

    identities = alignment["avg_identity"].tolist()
    conserved = diversity_summary["conserved_positions"].tolist()
    variable = diversity_summary["variable_positions"].tolist()
    diversity = diversity_summary["mean_diversity"].tolist()

    axes[0].bar(
        ['Type I', 'Type II'], identities,
        color=['#2E86AB', '#E76F51'], width=0.6
    )
    axes[0].set_ylabel('Mean Pairwise Identity')
    axes[0].set_ylim(0.5, 1.0)

    for i, v in enumerate(identities):
        axes[0].text(
            i, v + 0.01, f'{v:.3f}',
            ha='center', va='bottom', fontweight='bold'
        )

    axes[1].bar(
        [0, 1], conserved,
        width=0.4, color='#2E86AB',
        label='Conserved', align='center'
    )
    axes[1].bar(
        [0, 1], variable,
        width=0.4, bottom=conserved,
        color='#E76F51',
        label='Variable', align='center'
    )
    axes[1].set_xticks([0, 1])
    axes[1].set_xticklabels(['Type I', 'Type II'])
    axes[1].set_ylabel('Number of Positions')
    handles, labels = axes[1].get_legend_handles_labels()
    axes[1].legend(
        handles,
        labels,
        loc='upper left',
        bbox_to_anchor=(0.03, 0.98),
        ncol=2,
        frameon=False,
        fontsize=9,
        handlelength=1.5,
        columnspacing=1.5,
        borderaxespad=0.0
    )

    axes[2].bar(
        [0, 1], diversity,
        width=0.6, color=['#2E86AB', '#E76F51']
    )
    axes[2].set_xticks([0, 1])
    axes[2].set_xticklabels(['Type I', 'Type II'])
    axes[2].set_ylabel('Dominance-Complement Diversity')
    axes[2].set_ylim(0, 0.2)

    for i, v in enumerate(diversity):
        axes[2].text(
            i, v + 0.005, f'{v:.4f}',
            ha='center', va='bottom', fontweight='bold'
        )

    # Panel labels are deliberately placed inside each axes,
    # independent of the figure-level heading.
    for ax, label in zip(axes, ['A', 'B', 'C']):
        ax.text(
            0.0, 1.01, label,
            transform=ax.transAxes,
            ha='left', va='bottom',
            fontsize=11, fontweight='bold'
        )

    # One clean heading spanning the complete three-panel figure.
    fig.text(
        0.5, 0.985,
        'PRRSV Sequence Diversity and Alignment Characteristics',
        ha='center', va='top',
        fontsize=14, fontweight='bold'
    )

    # Separate figure note; it must not compete with the x-axis labels.
    fig.text(
        0.5, 0.012,
        'Dominance-complement diversity is a VAXINTAIC-defined comparative metric',
        ha='center', va='bottom',
        fontsize=9, style='italic'
    )

    # Fixed publication margins; do not call tight_layout after this.
    fig.subplots_adjust(
        left=0.075, right=0.985,
        top=0.88, bottom=0.16,
        wspace=0.30
    )

    for ext in ['pdf', 'tiff']:
        fig.savefig(
            OUTPUT_ROOT / 'figures' / f'Figure2_diversity.{ext}',
            dpi=300,
            bbox_inches=None
        )

    print('✓ Figure 2 saved')
    plt.close(fig)

def figure_3_candidate_landscape(screening, ranking):
    set_style()

    # Publication layout with a dedicated figure heading and
    # sufficient separation between axes and the explanatory note.
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.2))

    type_i = screening[screening['prrsv_type'] == 'Type I']
    type_ii = screening[screening['prrsv_type'] == 'Type II']

    axes[0].hist(
        type_i['epitope_score'], bins=50, alpha=0.5,
        color='#2E86AB', label=f'Type I (n={len(type_i)})'
    )
    axes[0].hist(
        type_ii['epitope_score'], bins=50, alpha=0.5,
        color='#E76F51', label=f'Type II (n={len(type_ii)})'
    )
    axes[0].set_xlabel('Peptide Score')
    axes[0].set_ylabel('Frequency')
    axes[0].legend(loc='upper right')

    type_i_rank = ranking[ranking['prrsv_type'] == 'Type I']
    type_ii_rank = ranking[ranking['prrsv_type'] == 'Type II']

    # One contiguous coordinate system for the complete screening landscape.
    screening_plot = screening.reset_index(drop=True).copy()
    screening_plot['plot_index'] = screening_plot.index

    type_i_plot = screening_plot[
        screening_plot['prrsv_type'] == 'Type I'
    ]
    type_ii_plot = screening_plot[
        screening_plot['prrsv_type'] == 'Type II'
    ]

    axes[1].scatter(
        type_i_plot['plot_index'],
        type_i_plot['epitope_score'],
        color='#2E86AB', alpha=0.3, s=5
    )
    axes[1].scatter(
        type_ii_plot['plot_index'],
        type_ii_plot['epitope_score'],
        color='#E76F51', alpha=0.3, s=5
    )

    screening_plot['key'] = (
        screening_plot['sequence'].astype(str)
        + '_'
        + screening_plot['start'].astype(str)
    )
    idx_map = dict(zip(screening_plot['key'], screening_plot.index))

    i_indices = []
    i_scores = []
    for _, r in type_i_rank.iterrows():
        key = str(r['sequence']) + '_' + str(r['start'])
        if key in idx_map:
            idx = idx_map[key]
            i_indices.append(idx)
            i_scores.append(screening_plot.loc[idx, 'epitope_score'])

    ii_indices = []
    ii_scores = []
    for _, r in type_ii_rank.iterrows():
        key = str(r['sequence']) + '_' + str(r['start'])
        if key in idx_map:
            idx = idx_map[key]
            ii_indices.append(idx)
            ii_scores.append(screening_plot.loc[idx, 'epitope_score'])

    if i_indices:
        axes[1].scatter(
            i_indices, i_scores,
            color='#2E86AB', edgecolor='black',
            s=50, marker='o', label='Top 15 Type I'
        )

    if ii_indices:
        axes[1].scatter(
            ii_indices, ii_scores,
            color='#E76F51', edgecolor='black',
            s=50, marker='^', label='Top 15 Type II'
        )

    axes[1].set_xlabel('Candidate Window Index')
    axes[1].set_ylabel('Peptide Score')
    axes[1].legend(
        loc='upper right',
        fontsize=8,
        frameon=False,
        handlelength=1.4,
        borderaxespad=0.4
    )

    axes[2].scatter(
        type_i_rank['start'],
        [1] * len(type_i_rank),
        color='#2E86AB', s=80, marker='o',
        label='Type I'
    )
    axes[2].scatter(
        type_ii_rank['start'],
        [0] * len(type_ii_rank),
        color='#E76F51', s=80, marker='^',
        label='Type II'
    )
    axes[2].set_xlabel('Genome Position')
    axes[2].set_yticks([0, 1])
    axes[2].set_yticklabels(['Type II', 'Type I'])
    axes[2].set_xlim(0, 600)
    axes[2].legend(
        loc='upper left',
        fontsize=8,
        frameon=False,
        handlelength=1.4,
        borderaxespad=0.4
    )

    # Clean, single figure-level heading.
    fig.text(
        0.5, 0.965,
        'Sequence-Derived Peptide Candidate Landscape and Selection',
        ha='center', va='top',
        fontsize=13, fontweight='bold'
    )

    # Independent panel labels; no set_title/suptitle interaction.
    for ax, label in zip(axes, ['A', 'B', 'C']):
        ax.text(
            0.0, 1.01, label,
            transform=ax.transAxes,
            ha='left', va='bottom',
            fontsize=11, fontweight='bold'
        )

    # Separate explanatory note from both x-axis labels.
    fig.text(
        0.5, 0.012,
        'Peptide scores combine hydrophilicity, flexibility, and conservation features',
        ha='center', va='bottom',
        fontsize=9, style='italic'
    )

    # Fixed margins; avoid tight_layout/bbox_inches='tight'
    # because they can compress the title/note into the panels.
    fig.subplots_adjust(
        left=0.075, right=0.985,
        top=0.88, bottom=0.16,
        wspace=0.30
    )

    for ext in ['pdf', 'tiff']:
        fig.savefig(
            OUTPUT_ROOT / 'figures' / f'Figure3_candidate_landscape.{ext}',
            dpi=300,
            bbox_inches=None
        )

    print('✓ Figure 3 saved')
    plt.close(fig)

def figure_4_construct_prioritization():
    set_style()
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    type_order = ['Type I', 'Type II']

    # Load authoritative VAXINTAIC outputs.
    construct_metadata = pd.read_csv(FILES['constructs'])
    cai_data = pd.read_csv(FILES['cai'])
    immunogenicity_data = pd.read_csv(FILES['immunogenicity'])
    nanoparticle_data = pd.read_csv(FILES['nanoparticle'])
    model_data = pd.read_csv(FILES['model'])

    # Panel A: Architecture diagram
    ax = axes[0, 0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis('off')
    comps = [
        (2.0, 0.8, 6.0, 2.0, 'Coding Sequence\n681 nt, 226 aa + stop', '#F4A261'),
        (0.5, 0.8, 1.5, 2.0, '5\' Region', '#2E86AB'),
        (8.0, 0.8, 1.5, 2.0, '3\' Region', '#2E86AB'),
    ]
    for x, y, w, h, label, color in comps:
        rect = FancyBboxPatch(
            (x, y), w, h, boxstyle='round,pad=0.1',
            facecolor=color, edgecolor='black', linewidth=1.5, alpha=0.7
        )
        ax.add_patch(rect)
        ax.text(
            x + w/2, y + h/2, label,
            ha='center', va='center', fontsize=9, fontweight='bold'
        )
    ax.text(
        5, 0.3,
        'Selected peptide candidates (15 per genotype) distributed along CDS',
        ha='center', va='center', fontsize=8, style='italic'
    )
    ax.set_title('A', loc='left', fontweight='bold')

    # Panel B: GC/CAI
    ax = axes[0, 1]
    x = [0, 1]

    gc = [
        float(construct_metadata.loc[
            construct_metadata['prrsv_type'] == t, 'GC_percent'
        ].iloc[0])
        for t in type_order
    ]
    cai = [
        float(cai_data.loc[
            cai_data['construct_id'] ==
            construct_metadata.loc[
                construct_metadata['prrsv_type'] == t, 'construct_id'
            ].iloc[0],
            'CAI'
        ].iloc[0])
        for t in type_order
    ]

    w = 0.35
    ax.bar(
        [p - w/2 for p in x], gc, w,
        color=['#2E86AB', '#E76F51'], label='GC%'
    )
    ax.bar(
        [p + w/2 for p in x], [c * 50 for c in cai], w,
        color=['#2E86AB', '#E76F51'], alpha=0.5,
        label='CAI (×50)', hatch='//'
    )
    ax.set_xticks(x)
    ax.set_xticklabels(type_order)
    ax.set_ylabel('Percentage')
    ax.set_title('B', loc='left', fontweight='bold')
    ax.legend(loc='upper right', fontsize=8)

    for i, v in enumerate(gc):
        ax.text(
            i - w/2, v + 0.5, f'{v:.2f}%',
            ha='center', va='bottom', fontsize=8
        )
    for i, v in enumerate(cai):
        ax.text(
            i + w/2, v * 50 + 0.5, f'{v:.3f}',
            ha='center', va='bottom', fontsize=8
        )

    # Panel C: Immunogenicity + in-silico delivery estimates
    ax = axes[1, 0]
    metrics = [
        'Immunogenicity',
        'Encapsulation\nEstimate',
        'Delivery\nStability'
    ]

    t1 = [
        float(immunogenicity_data.loc[
            immunogenicity_data['prrsv_type'] == 'Type I',
            'immunogenicity'
        ].mean()),
        float(nanoparticle_data.loc[
            nanoparticle_data['prrsv_type'] == 'Type I',
            'encapsulation_efficiency'
        ].iloc[0]),
        float(nanoparticle_data.loc[
            nanoparticle_data['prrsv_type'] == 'Type I',
            'delivery_stability'
        ].iloc[0])
    ]

    t2 = [
        float(immunogenicity_data.loc[
            immunogenicity_data['prrsv_type'] == 'Type II',
            'immunogenicity'
        ].mean()),
        float(nanoparticle_data.loc[
            nanoparticle_data['prrsv_type'] == 'Type II',
            'encapsulation_efficiency'
        ].iloc[0]),
        float(nanoparticle_data.loc[
            nanoparticle_data['prrsv_type'] == 'Type II',
            'delivery_stability'
        ].iloc[0])
    ]

    x = [0, 1, 2]
    ax.bar(
        [p - 0.2 for p in x], t1, 0.4,
        color='#2E86AB', label='Type I'
    )
    ax.bar(
        [p + 0.2 for p in x], t2, 0.4,
        color='#E76F51', label='Type II'
    )
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylabel('Score / Estimate')
    ax.set_title('C', loc='left', fontweight='bold')
    ax.legend(loc='upper right', fontsize=8)
    ax.set_ylim(0, 1.0)

    # Panel D: Integrated computational prioritization index
    ax = axes[1, 1]
    x = [0, 1]
    vi = [
        float(model_data.loc[
            model_data['prrsv_type'] == t, 'vaccine_index'
        ].iloc[0])
        for t in type_order
    ]

    bars = ax.bar(
        x, vi, width=0.6,
        color=['#2E86AB', '#E76F51']
    )
    ax.set_xticks(x)
    ax.set_xticklabels(type_order)
    ax.set_ylabel('Vaccine Prioritization Index')
    ax.set_title('D', loc='left', fontweight='bold')
    ax.set_ylim(0, 1.0)

    for bar, v in zip(bars, vi):
        ax.text(
            bar.get_x() + bar.get_width()/2,
            v + 0.01,
            f'{v:.3f}',
            ha='center', va='bottom', fontweight='bold'
        )

    fig.suptitle(
        'mRNA Construct Architecture and Integrated Computational Prioritization',
        fontsize=14, fontweight='bold'
    )
    fig.text(
        0.5, 0.01,
        'All scores are VAXINTAIC computational prioritization measures, '
        'not experimentally demonstrated efficacy',
        ha='center', va='center', fontsize=9, style='italic'
    )
    plt.tight_layout()
    plt.subplots_adjust(top=0.93, bottom=0.06)

    for ext in ['pdf', 'tiff']:
        plt.savefig(
            OUTPUT_ROOT / 'figures' /
            f'Figure4_construct_prioritization.{ext}',
            dpi=300, bbox_inches='tight'
        )

    print('✓ Figure 4 saved')
    plt.close()

# ============================================================================
# GENERATE TABLES
# ============================================================================

def generate_tables(screening, ranking):
    (OUTPUT_ROOT / 'tables').mkdir(parents=True, exist_ok=True)

    # Table 1: Dataset, alignment, and diversity
    alignment = pd.read_csv(FILES['alignment'])
    diversity_summary = pd.read_csv(FILES['diversity'])

    alignment_by_type = alignment.set_index('type')
    diversity_by_type = diversity_summary.set_index('prrsv_type')

    rows1 = [
        ['Number of sequences (n)',
         str(int(alignment_by_type.loc['Type I', 'num_sequences'])),
         str(int(alignment_by_type.loc['Type II', 'num_sequences']))],
        ['Sequence length (nt)',
         '600–606',
         '600–606'],
        ['Alignment length (nt)',
         str(int(alignment_by_type.loc['Type I', 'alignment_length'])),
         str(int(alignment_by_type.loc['Type II', 'alignment_length']))],
        ['Mean pairwise identity',
         f"{float(alignment_by_type.loc['Type I', 'avg_identity']):.4f}",
         f"{float(alignment_by_type.loc['Type II', 'avg_identity']):.4f}"],
        ['Mean gap proportion',
         f"{float(alignment_by_type.loc['Type I', 'gap_fraction']):.4f}",
         f"{float(alignment_by_type.loc['Type II', 'gap_fraction']):.4f}"],
        ['Conserved positions',
         str(int(diversity_by_type.loc['Type I', 'conserved_positions'])),
         str(int(diversity_by_type.loc['Type II', 'conserved_positions']))],
        ['Variable positions',
         str(int(diversity_by_type.loc['Type I', 'variable_positions'])),
         str(int(diversity_by_type.loc['Type II', 'variable_positions']))],
        ['Dominance-complement diversity',
         f"{float(diversity_by_type.loc['Type I', 'mean_diversity']):.4f}",
         f"{float(diversity_by_type.loc['Type II', 'mean_diversity']):.4f}"],
    ]

    df1 = pd.DataFrame(rows1, columns=['Parameter', 'Type I', 'Type II'])
    df1.to_csv(OUTPUT_ROOT / 'tables' / 'Table1_dataset_alignment.csv', index=False)
    print('✓ Table 1 saved')

    # Table 2: Top five unique candidate peptide sequences per genotype.
    # The complete 30-record selected candidate set is retained in
    # Supplementary_Table_S2_top_candidates.csv.
    t1 = ranking[ranking['prrsv_type'] == 'Type I'].drop_duplicates(
        subset=['sequence'], keep='first'
    ).head(5)

    t2 = ranking[ranking['prrsv_type'] == 'Type II'].drop_duplicates(
        subset=['sequence'], keep='first'
    ).head(5)

    df2 = pd.concat([t1, t2], ignore_index=True)
    df2 = df2.rename(columns={
        'prrsv_type': 'Genotype',
        'epitope_score': 'Peptide Score'
    })

    cols = [
        'Genotype', 'sequence', 'start', 'conservation',
        'hydrophilicity', 'flexibility', 'Peptide Score'
    ]
    cols = [c for c in cols if c in df2.columns]
    df2 = df2[cols]

    for col in ['conservation', 'hydrophilicity', 'flexibility', 'Peptide Score']:
        if col in df2.columns:
            df2[col] = df2[col].map(lambda x: f'{float(x):.3f}')

    df2.to_csv(
        OUTPUT_ROOT / 'tables' / 'Table2_peptide_candidates.csv',
        index=False
    )
    print('✓ Table 2 saved')

    # Table 3: Construct and computational prioritization
    construct_metadata = pd.read_csv(FILES['constructs'])
    cai_data = pd.read_csv(FILES['cai'])
    immunogenicity_data = pd.read_csv(FILES['immunogenicity'])
    nanoparticle_data = pd.read_csv(FILES['nanoparticle'])
    model_data = pd.read_csv(FILES['model'])

    type_order = ['Type I', 'Type II']

    rows3 = [
        ['Genotype', 'Type I', 'Type II'],
        ['Nucleotide length (nt)'] + [
            str(int(construct_metadata.loc[
                construct_metadata['prrsv_type'] == t, 'length_nt'
            ].iloc[0]))
            for t in type_order
        ],
        ['Amino acid length (aa) + terminal stop', '226 + stop', '226 + stop'],
        ['GC content (%)'] + [
            f"{float(construct_metadata.loc[
                construct_metadata['prrsv_type'] == t, 'GC_percent'
            ].iloc[0]):.2f}"
            for t in type_order
        ],
        ['Codon Adaptation Index (CAI)'] + [
            f"{float(cai_data.loc[
                cai_data['construct_id'] ==
                construct_metadata.loc[
                    construct_metadata['prrsv_type'] == t, 'construct_id'
                ].iloc[0],
                'CAI'
            ].iloc[0]):.3f}"
            for t in type_order
        ],
        ['Mean immunogenicity score'] + [
            f"{float(immunogenicity_data.loc[
                immunogenicity_data['prrsv_type'] == t,
                'immunogenicity'
            ].mean()):.3f}"
            for t in type_order
        ],
        ['Encapsulation-efficiency estimate'] + [
            f"{float(nanoparticle_data.loc[
                nanoparticle_data['prrsv_type'] == t,
                'encapsulation_efficiency'
            ].iloc[0]):.3f}"
            for t in type_order
        ],
        ['Delivery-stability estimate'] + [
            f"{float(nanoparticle_data.loc[
                nanoparticle_data['prrsv_type'] == t,
                'delivery_stability'
            ].iloc[0]):.3f}"
            for t in type_order
        ],
        ['Delivery index'] + [
            f"{float(nanoparticle_data.loc[
                nanoparticle_data['prrsv_type'] == t,
                'delivery_index'
            ].iloc[0]):.3f}"
            for t in type_order
        ],
        ['Model score'] + [
            f"{float(model_data.loc[
                model_data['prrsv_type'] == t,
                'predicted_model_score'
            ].iloc[0]):.3f}"
            for t in type_order
        ],
        ['Vaccine prioritization index'] + [
            f"{float(model_data.loc[
                model_data['prrsv_type'] == t,
                'vaccine_index'
            ].iloc[0]):.3f}"
            for t in type_order
        ],
    ]

    df3 = pd.DataFrame(
        rows3,
        columns=['Parameter', 'Type I Construct', 'Type II Construct']
    )
    df3.to_csv(
        OUTPUT_ROOT / 'tables' / 'Table3_construct_prioritization.csv',
        index=False
    )
    print('✓ Table 3 saved')

# ============================================================================
# GENERATE SUPPLEMENTARY
# ============================================================================

def generate_supplementary(screening, ranking, constructs):
    sup = OUTPUT_ROOT / 'supplementary'
    sup.mkdir(parents=True, exist_ok=True)

    screening.to_csv(sup / 'Supplementary_Table_S1_full_screening.csv', index=False)
    ranking.to_csv(sup / 'Supplementary_Table_S2_top_candidates.csv', index=False)
    constructs.to_csv(sup / 'Supplementary_Table_S3_construct_metadata.csv', index=False)
    print('✓ Supplementary tables saved')

# ============================================================================
# MAIN
# ============================================================================

def main():
    print('\n' + '=' * 70)
    print('VAXINTAIC TVJ Publication Generator (Production Version)')
    print('=' * 70 + '\n')

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    (OUTPUT_ROOT / 'figures').mkdir(exist_ok=True)
    (OUTPUT_ROOT / 'tables').mkdir(exist_ok=True)
    (OUTPUT_ROOT / 'supplementary').mkdir(exist_ok=True)

    try:
        data, hashes = load_data()
    except (FileNotFoundError, ValueError) as e:
        print('\n' + '!' * 70)
        print('ERROR: Publication generation FAILED')
        print('!' * 70)
        print(f'\n{str(e)}')
        print('\nThe generator has stopped. No synthetic data has been created.')
        sys.exit(1)

    verify_data(data)

    print('\nGenerating figures...')
    figure_1_workflow()
    figure_2_diversity()
    figure_3_candidate_landscape(data['screening'], data['ranking'])
    figure_4_construct_prioritization()

    print('\nGenerating tables...')
    generate_tables(data['screening'], data['ranking'])

    print('\nGenerating supplementary...')
    generate_supplementary(data['screening'], data['ranking'], data['constructs'])

    # Validation report
    report_path = OUTPUT_ROOT / 'VALIDATION_REPORT.txt'
    with open(report_path, 'w') as f:
        f.write('VAXINTAIC PUBLICATION DATA VALIDATION REPORT\n')
        f.write('=' * 70 + '\n\n')
        f.write(f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
        f.write('STATUS: PUBLICATION DATA VERIFIED\n')
        f.write('No synthetic/fabricated data used.\n')
    print('✓ Validation report saved')

    print('\n' + '=' * 70)
    print('✅ PUBLICATION PACKAGE GENERATED SUCCESSFULLY')
    print('=' * 70)
    print(f'\n📁 Output: {OUTPUT_ROOT}')

if __name__ == '__main__':
    main()
