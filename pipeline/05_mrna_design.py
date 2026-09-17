import os

# SCIENTIFIC VALIDITY CONFIGURATION
#
# Host organism must be explicit because codon optimization is host-specific.
# Existing Stage 05 behavior indicates Sus scrofa; Stage 11 must be reconciled
# against this setting before final scientific claims are made.
TARGET_HOST_ORGANISM = os.environ.get(
    "VAXINTAIC_TARGET_HOST_ORGANISM",
    "Sus scrofa",
)

#!/usr/bin/env python3
"""
VAXINTAIC – Stage 05 : Current Public mRNA Construct Design

Builds codon-optimized mRNA constructs from top-ranked epitopes (Stage 04).
Automatically translates nucleotide consensus sequences into amino acids.
Outputs codon-optimized DNA (mRNA equivalent), with QC metrics.
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# ────────────────────────────────────────────────
# CONFIG
# ────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
TOP_EPITOPES_FILE = BASE_DIR / "data" / "candidates" / "top_epitopes_public_current.csv"
OUT_DIR = BASE_DIR / "data" / "mrna"
os.makedirs(OUT_DIR, exist_ok=True)

# Deterministic local RNG for reproducible codon selection.
RNG_SEED = int(os.environ.get("VAXINTAIC_STAGE05_SEED", "42"))
RNG = np.random.default_rng(RNG_SEED)

if not os.path.exists(TOP_EPITOPES_FILE):
    raise SystemExit("❌ Missing top_epitopes_public_current.csv — run Stage 04 first.")

# Host-specific codon usage is loaded from the validated
# Sus scrofa reference in config.codon_reference.
from config.codon_reference import (
    TARGET_HOST_ORGANISM,
    SUS_SCROFA_CODONS_BY_AA,
    SUS_SCROFA_CODON_USAGE,
    CODON_REFERENCE_STATUS,
    CODON_REFERENCE_SOURCE,
)

if TARGET_HOST_ORGANISM != "Sus scrofa":
    raise SystemExit("Target host must be Sus scrofa.")

if not SUS_SCROFA_CODONS_BY_AA:
    raise SystemExit(
        "mRNA design BLOCKED: documented Sus scrofa codon reference unavailable."
    )

CODON_USAGE = SUS_SCROFA_CODONS_BY_AA

# ────────────────────────────────────────────────
# Helper functions
# ────────────────────────────────────────────────
def codon_optimize(peptide):
    """Convert amino acids to DNA using Sus scrofa relative-adaptiveness weights."""
    dna_seq = ""
    for aa in peptide:
        if aa not in CODON_USAGE:
            raise ValueError(
                f"Unsupported amino-acid residue {aa}; refusing to generate NNN."
            )

        codons = CODON_USAGE[aa]
        weights = np.array(
            [SUS_SCROFA_CODON_USAGE[codon] for codon in codons],
            dtype=float,
        )

        if not np.isfinite(weights).all() or weights.sum() <= 0:
            raise ValueError(
                f"Invalid Sus scrofa codon weights for amino acid {aa}."
            )

        weights /= weights.sum()
        dna_seq += RNG.choice(codons, p=weights)

    return dna_seq

def calc_gc(seq):
    return round((seq.count('G') + seq.count('C')) / len(seq) * 100, 2)

def calc_cai(seq):
    """Calculate CAI as the geometric mean of validated Sus scrofa weights."""
    weights = []

    for i in range(0, len(seq) - 3, 3):
        codon = seq[i:i+3].upper().replace("T", "U")

        # Exclude terminal stop codon from CAI.
        if codon in {"UAA", "UAG", "UGA"}:
            continue

        weight = SUS_SCROFA_CODON_USAGE.get(codon)
        if weight is None or not (0 < float(weight) <= 1):
            raise ValueError(
                f"Invalid or missing Sus scrofa CAI weight for codon {codon}."
            )

        weights.append(float(weight))

    if not weights:
        return 0.0

    return round(float(np.exp(np.mean(np.log(weights)))), 3)

# ────────────────────────────────────────────────
# LOAD EPITOPES
# ────────────────────────────────────────────────
epitopes = pd.read_csv(TOP_EPITOPES_FILE)
print(f"📥 Loaded {len(epitopes)} top epitopes from: {TOP_EPITOPES_FILE}")

# ────────────────────────────────────────────────
# BUILD CONSTRUCTS PER TYPE
# ────────────────────────────────────────────────
construct_records = []
metadata_rows = []

for prrsv_type, group in epitopes.groupby("prrsv_type"):
    joined_peptides = ''.join(group["consensus_seq"].tolist())

    # 🧠 Auto-detect if sequence looks like DNA (A/T/G/C/U)
    if set(joined_peptides.upper()) <= set("ATGCU"):
        try:
            joined_peptides = str(Seq(joined_peptides).translate(to_stop=True))
            print(f"🔄 Translated {prrsv_type} consensus to amino acids ({len(joined_peptides)} aa).")
        except Exception as e:
            print(f"⚠️ Translation failed for {prrsv_type}: {e}")

    # Codon optimization
    dna_seq = codon_optimize(joined_peptides)
    dna_seq = "ATG" + dna_seq + "TAA"  # add start/stop codons

    # QC metrics
    gc_content = calc_gc(dna_seq)
    cai = calc_cai(dna_seq)
    seq_len = len(dna_seq)
    construct_id = f"{prrsv_type.replace(' ', '_')}_Construct_1"

    construct_records.append(SeqRecord(
        Seq(dna_seq),
        id=construct_id,
        description=f"{prrsv_type} optimized mRNA construct"
    ))

    metadata_rows.append({
        "construct_id": construct_id,
        "prrsv_type": prrsv_type,
        "num_epitopes": len(group),
        "length_nt": seq_len,
        "GC_percent": gc_content,
        "CAI": cai
    })

    # Print summary per construct
    print(f"🧬 {construct_id}: length={seq_len}nt, GC%={gc_content}, CAI={cai}")

# ────────────────────────────────────────────────
# SAVE OUTPUTS
# ────────────────────────────────────────────────
fasta_out = os.path.join(OUT_DIR, "mrna_constructs_public_current.fasta")
meta_out  = os.path.join(OUT_DIR, "mrna_construct_metadata_public_current.csv")

SeqIO.write(construct_records, fasta_out, "fasta")
pd.DataFrame(metadata_rows).to_csv(meta_out, index=False)

print(f"\n✅ mRNA constructs saved → {fasta_out}")
print(f"📄 Construct metadata → {meta_out}")
print("🧬 Stage 05 completed successfully.")

