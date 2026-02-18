#!/usr/bin/env python3

import argparse
from pathlib import Path
from Bio import SeqIO
import pandas as pd

# -----------------------------
# CLI
# -----------------------------
parser = argparse.ArgumentParser(
    description="Verify expected bases at specific gene-relative positions"
)
parser.add_argument(
    "--fasta",
    required=True,
    help="FASTA file to verify (e.g. ATCC19977_WT.fasta or a mutant FASTA)"
)
parser.add_argument(
    "--outdir",
    required=False,
    default=".",
    help="Output directory for verification results (default: current directory)"
)
args = parser.parse_args()

FASTA_PATH = Path(args.fasta).resolve()
OUTDIR = Path(args.outdir).resolve()
OUTDIR.mkdir(parents=True, exist_ok=True)

fasta_name = FASTA_PATH.name
if fasta_name.endswith(".fasta"):
    fasta_name = fasta_name[:-6]  # strip ".fasta"

OUTFILE = OUTDIR / f"{fasta_name}_site_verification.tsv"

# -----------------------------
# Expected sites (1-based)
# These are gene-relative expectations
# -----------------------------
EXPECTED = {
    "rrl": {
        "start": 1464208,   # 1-based
        "sites": {
            2269: "A",
            2270: "A",
            2271: "A",
            2281: "G",
            2293: "A",
        },
    },
    "rrs": {
        "start": 1462398,   # 1-based
        "sites": {
            1373: "T",
            1375: "A",
            1376: "C",
            1458: "G",
        },
    },
    "erm41": {
        "start": 2345955,   # 1-based
        "sites": {
            19: "C",
            28: "T",
        },
    },
}

# -----------------------------
# Load FASTA
# -----------------------------
if not FASTA_PATH.exists():
    raise FileNotFoundError(f"FASTA not found: {FASTA_PATH}")

record = SeqIO.read(FASTA_PATH, "fasta")
seq = str(record.seq).upper()
seq_id = record.id
seq_len = len(seq)

# -----------------------------
# Main verification
# -----------------------------
rows = []

for gene, info in EXPECTED.items():
    gene_start = info["start"]

    for gene_pos, expected_base in info["sites"].items():
        ref_pos = gene_start + gene_pos - 1  # convert to reference coordinate

        if ref_pos < 1 or ref_pos > seq_len:
            observed = "OUT_OF_RANGE"
            match = False
        else:
            observed = seq[ref_pos - 1]  # 1-based → 0-based
            match = (observed == expected_base)

        rows.append({
            "fasta": FASTA_PATH.name,
            "sequence_id": seq_id,
            "gene": gene,
            "gene_position_1based": gene_pos,
            "ref_position_1based": ref_pos,
            "expected_base": expected_base,
            "observed_base": observed,
            "match": match,
        })

df = pd.DataFrame(rows)
df.to_csv(OUTFILE, sep="\t", index=False)

print(f"\nSite verification written to: {OUTFILE}\n")
print(df.to_string(index=False))

