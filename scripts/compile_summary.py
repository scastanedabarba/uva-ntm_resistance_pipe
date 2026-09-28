#!/usr/bin/env python3

import argparse
import os
import statistics
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

# -----------------------------------
# CLI
# -----------------------------------
parser = argparse.ArgumentParser(
    description="Step3 compile: BLAST top hits + complete target-site evidence across isolates + combined status tables."
)
parser.add_argument(
    "-o", "--outdir", required=True,
    help="Pipeline output root (contains per-isolate subdirs + summary/)"
)
parser.add_argument(
    "-l", "--linelist", required=True,
    help="Linelist (CSV or TSV). Isolate must be in first column. Run column is ignored."
)
parser.add_argument(
    "-f", "--targets-fasta", required=True,
    help="Targets FASTA used in Step 1 BLAST (e.g., scripts/nucleotide.fna)"
)
parser.add_argument(
    "-r", "--ref-contig", default="CU458896.1",
    help="Reference contig name (only used for depth lookup)"
)
args = parser.parse_args()

OUTDIR = os.path.abspath(args.outdir)
LINELIST = os.path.abspath(args.linelist)
TARGETS_FASTA = os.path.abspath(args.targets_fasta)
REF_CONTIG = args.ref_contig.strip()

os.makedirs(os.path.join(OUTDIR, "summary"), exist_ok=True)
os.makedirs(os.path.join(OUTDIR, "status"), exist_ok=True)

# -----------------------------------
# Target gene coordinates (MATCH 02_map_call.slurm)
# 1-based inclusive
# -----------------------------------
RRS_START, RRS_END = 1462398, 1463901
RRL_START, RRL_END = 1464208, 1467319
ERM41_START, ERM41_END = 2345955, 2346476

GENE_STARTS = {
    "rrs": RRS_START,
    "rrl": RRL_START,
    "erm41": ERM41_START,
}

# Sites-of-interest (gene coords, 1-based)
SITES: List[Tuple[str, int]] = [
    ("rrl", 2269), ("rrl", 2270), ("rrl", 2271), ("rrl", 2281), ("rrl", 2293),
    ("erm41", 19), ("erm41", 28),
    ("rrs", 1373), ("rrs", 1375), ("rrs", 1376), ("rrs", 1458),
]

# Reference alleles at predefined sites of interest.
# Used to populate complete site evidence even when no variant is called.
REFERENCE_BASES: Dict[Tuple[str, int], str] = {
    ("rrl", 2269): "A",
    ("rrl", 2270): "A",
    ("rrl", 2271): "A",
    ("rrl", 2281): "G",
    ("rrl", 2293): "A",
    ("erm41", 19): "C",
    ("erm41", 28): "T",
    ("rrs", 1373): "T",
    ("rrs", 1375): "A",
    ("rrs", 1376): "C",
    ("rrs", 1458): "G",
}

# -----------------------------------
# erm41 truncation metrics (coverage-based; Criterion A prep)
# -----------------------------------
ERM41_LEFT_FLANK = (20, 140)     # gene positions (1-based, inclusive)
ERM41_RIGHT_FLANK = (450, 560)   # gene positions (1-based, inclusive)
ERM41_DEL_RANGE = (159, 432)     # gene positions (1-based, inclusive)
ERM41_CALLABLE_MIN_FLANK_MED = 20  # "callable" threshold on flank median depth

# -----------------------------------
# Helpers
# -----------------------------------
def read_isolates_from_linelist(path: str) -> List[str]:
    """
    Accept CSV or TSV. Header optional.
    Uses first column as isolate.
    """
    isolates: List[str] = []
    with open(path, "r") as f:
        for line in f:
