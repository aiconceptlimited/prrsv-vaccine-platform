#!/usr/bin/env python3
from pathlib import Path
"""
VAXINTAIC v3 – Stage 02 (Enhanced):
Multiple Sequence Alignment + Pairwise Identity Analysis

Aligns PRRSV Type I and Type II sequences using MAFFT and computes:
 - Number of sequences
 - Alignment length
 - Gap fraction
 - Mean pairwise nucleotide identity
"""

import os, subprocess, datetime, pandas as pd
from Bio import AlignIO

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "sequences"
ALIGN_DIR = BASE_DIR / "data" / "alignments"
META_FILE = os.path.join(DATA_DIR, "metadata_v3.csv")
FASTA_FILE = os.path.join(DATA_DIR, "prrsv_orf5or6_clean.fasta")
SUMMARY_FILE = os.path.join(ALIGN_DIR, "alignment_summary_v3.csv")

os.makedirs(ALIGN_DIR, exist_ok=True)

print(f"[{datetime.datetime.now()}] 🚀 Starting PRRSV alignment + identity analysis stage…")

# ────────────────────────────────────────────────
# Load metadata and group by type
# ────────────────────────────────────────────────
meta = pd.read_csv(META_FILE)
types = meta["prrsv_type"].unique().tolist()
print(f"🧩 Sequence distribution by PRRSV type: {dict(meta['prrsv_type'].value_counts())}")

# ────────────────────────────────────────────────
# Helper: run MAFFT
# ────────────────────────────────────────────────
def run_mafft(input_fasta, output_fasta):
    cmd = ["mafft", "--auto", "--thread", "4", input_fasta]
    subprocess.run(cmd, stdout=open(output_fasta, "w"), stderr=subprocess.DEVNULL, check=True)

# ────────────────────────────────────────────────
# Helper: extract sequences by type
# ────────────────────────────────────────────────
def extract_type_fasta(prrsv_type):
    from Bio import SeqIO
    tmp_path = os.path.join(ALIGN_DIR, f"temp_{prrsv_type}.fasta")
    ids = meta.loc[meta["prrsv_type"] == prrsv_type, "id"].tolist()
    records = [r for r in SeqIO.parse(FASTA_FILE, "fasta") if r.id in ids]
    if not records:
        return None
    SeqIO.write(records, tmp_path, "fasta")
    return tmp_path

# ────────────────────────────────────────────────
# Helper: compute mean pairwise identity
# ────────────────────────────────────────────────
def compute_mean_identity(alignment):
    """Compute mean pairwise identity (ignoring gaps)."""
    from itertools import combinations
    import numpy as np

    num_seq = len(alignment)
    if num_seq < 2:
        return 1.0

    seqs = [str(rec.seq) for rec in alignment]
    total_identity, total_pairs = 0, 0

    for s1, s2 in combinations(seqs, 2):
        matches, valid = 0, 0
        for a, b in zip(s1, s2):
            if a == "-" or b == "-":
                continue
            valid += 1
            if a == b:
                matches += 1
        if valid > 0:
            total_identity += matches / valid
            total_pairs += 1

    return round(total_identity / total_pairs, 4) if total_pairs else 1.0

# ────────────────────────────────────────────────
# Align + analyze each genotype
# ────────────────────────────────────────────────
summary_rows = []
for t in ["Type I", "Type II"]:
    tmp = extract_type_fasta(t)
    if tmp is None:
        print(f"⚠️ No sequences found for {t}. Skipping.")
        continue

    out = os.path.join(ALIGN_DIR, f"prrsv_{t.replace(' ', '').lower()}_aligned.fasta")
    print(f"🔬 Aligning {t} sequences ({len(meta[meta['prrsv_type']==t])} seqs)…")

    try:
        run_mafft(tmp, out)
        aln = AlignIO.read(out, "fasta")
        aln_len = aln.get_alignment_length()
        num_seq = len(aln)
        gap_fraction = sum(rec.seq.count('-') for rec in aln) / (aln_len * num_seq)
        mean_identity = compute_mean_identity(aln)

        summary_rows.append({
            "type": t,
            "num_sequences": num_seq,
            "alignment_length": aln_len,
            "gap_fraction": round(gap_fraction, 4),
            "avg_identity": mean_identity
        })
        print(f"✅ {t} alignment complete → {out}")
        print(f"   ↳ Mean identity: {mean_identity:.3f}, Gaps: {gap_fraction:.3f}")
    except Exception as e:
        print(f"❌ Alignment failed for {t}: {e}")

# ────────────────────────────────────────────────
# Save summary
# ────────────────────────────────────────────────
if summary_rows:
    pd.DataFrame(summary_rows).to_csv(SUMMARY_FILE, index=False)
    print(f"📊 Alignment summary saved → {SUMMARY_FILE}")
else:
    print("⚠️ No alignments were performed.")

print("🧬 Stage 02 completed successfully.")

