#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# make_site_mutants.sh
#
# Creates ATCC19977 WT + per-site mutant FASTAs, including:
#
#   - rrl target-site mutants
#   - rrs target-site mutants
#   - erm41 C19T
#   - erm41 T28C
#   - literature-described truncated erm41 allele:
#       Δ64–65
#       Δ159–432
#
# Also creates additional erm-gene simulation controls:
#
#   - ATCC35855_erm39
#       Unmodified ATCC35855 genome containing native erm39
#
#   - CCUG47445_erm55_plasmid
#   - CCUG47445_erm55_transposon
#   - CCUG47445_erm55_chromosome
#
#     Each is based on the complete CCUG47445 chromosome
#     (CP007220.1) with one of the exact erm55 query sequences
#     from references/nucleotide.fna inserted at the same
#     synthetic location.
#
# The erm41 truncation simulation includes BOTH deletions, while
# downstream interpretation uses the large Δ159–432 deletion as
# the robust coverage-based truncation signal.
#
# ATCC19977 site verification is run only on ATCC19977-derived
# FASTAs because verify_sites.py uses ATCC19977 genomic coordinates.
#
# Outputs:
#
#   <workdir>/atcc_dataset/
#     refs/
#     vcfs/
#     tmp/
#       truth_table.tsv
#       truth_sv.tsv
#     site_verification/
#       <ATCC19977 fasta>_site_verification.tsv
#       site_verification_summary.tsv
# ------------------------------------------------------------

usage() {
  cat <<'EOF'
USAGE:
  bash make_site_mutants.sh --workdir <path> [options]

REQUIRED:
  --workdir <dir>         Working directory

OPTIONAL:
  --chrom <name>          ATCC19977 contig name (default: CU458896.1)
  --no_skip_same          Do NOT skip when ALT==REF
  -h, --help              Show help
EOF
}

WORKDIR=""
CHROM="CU458896.1"
NO_SKIP_SAME=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir)
      WORKDIR="${2:-}"
      shift 2
      ;;
    --chrom)
      CHROM="${2:-}"
      shift 2
      ;;
    --no_skip_same)
      NO_SKIP_SAME=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

[[ -z "$WORKDIR" ]] && {
  echo "ERROR: --workdir required" >&2
  exit 1
}

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

source "$REPO_ROOT/config/setup_environment.sh"

REF_DEFAULT="$REPO_ROOT/references/ATCC19977.fasta"
REF="${REF_FASTA:-$REF_DEFAULT}"

ATCC35855_REF="$REPO_ROOT/references/ATCC35855.fasta"
CCUG47445_REF="$REPO_ROOT/references/CCUG47445.fasta"
NUCLEOTIDE_REF="$REPO_ROOT/references/nucleotide.fna"

OUTDIR="${WORKDIR%/}/atcc_dataset"

[[ -f "$REF" ]] || {
  echo "ERROR: ATCC19977 reference FASTA not found: $REF" >&2
  exit 1
}

[[ -f "$ATCC35855_REF" ]] || {
  echo "ERROR: ATCC35855 reference FASTA not found: $ATCC35855_REF" >&2
  exit 1
}

[[ -f "$CCUG47445_REF" ]] || {
  echo "ERROR: CCUG47445 reference FASTA not found: $CCUG47445_REF" >&2
  exit 1
}

[[ -f "$NUCLEOTIDE_REF" ]] || {
  echo "ERROR: nucleotide reference FASTA not found: $NUCLEOTIDE_REF" >&2
  exit 1
}

# ------------------------------------------------------------
# Environment
# ------------------------------------------------------------

setup_mutant_env
setup_python_env

mkdir -p "$OUTDIR"/{refs,vcfs,tmp}

# ------------------------------------------------------------
# Gene starts in ATCC19977 reference
# 1-based genomic coordinates
# ------------------------------------------------------------

RRL_START=1464208
RRS_START=1462398
ERM41_START=2345955

# ------------------------------------------------------------
# Sites of interest
# Gene-relative, 1-based
# ------------------------------------------------------------

RRL_SITES=(2269 2270 2271 2281 2293)
RRS_SITES=(1373 1375 1376 1458)

# Convert a 1-based gene-relative coordinate to a
# 1-based reference-genome coordinate.
gene_to_refpos() {
  echo $(( $1 + $2 - 1 ))
}

# ------------------------------------------------------------
# Index ATCC19977 reference
# ------------------------------------------------------------

[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

if ! cut -f1 "${REF}.fai" | grep -Fxq "$CHROM"; then
  echo "ERROR: Contig '$CHROM' not found in FASTA index" >&2
  echo "Available contigs (first 20):" >&2
  cut -f1 "${REF}.fai" | head -n 20 >&2
  exit 1
fi

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

get_ref_base() {
  samtools faidx "$REF" "${CHROM}:$1-$1" \
    | awk 'NR==2 {print toupper($0)}'
}

flip_base() {
  case "$1" in
    A) echo G ;;
    G) echo A ;;
    C) echo T ;;
    T) echo C ;;
    *) echo A ;;
  esac
}

# ------------------------------------------------------------
# WT reference
# ------------------------------------------------------------

WT="$OUTDIR/refs/ATCC19977_WT.fasta"

cp -f "$REF" "$WT"
samtools faidx "$WT"

# ------------------------------------------------------------
# SNP truth table
# ------------------------------------------------------------

TRUTH="$OUTDIR/tmp/truth_table.tsv"

printf "Dataset\tCHROM\tPOS_ref\tREF\tALT\n" > "$TRUTH"

# ------------------------------------------------------------
# SNP mutant generation
# ------------------------------------------------------------

make_mutant() {
  local label="$1"
  local pos="$2"
  local forced_alt="${3:-}"

  local ref
  local alt

  ref="$(get_ref_base "$pos")"

  [[ -z "$ref" ]] && {
    echo "WARNING: no REF at $pos"
    return 0
  }

  if [[ -n "$forced_alt" ]]; then
    alt="$forced_alt"
  else
    alt="$(flip_base "$ref")"
  fi

  if [[ "$alt" == "$ref" && "$NO_SKIP_SAME" -eq 0 ]]; then
    echo \
      "WARNING: ALT==REF at ${CHROM}:${pos} for ${label}. Skipping." \
      >&2
    return 0
  fi

  local vcf="$OUTDIR/vcfs/${label}.vcf"
  local vcfgz="${vcf}.gz"

  {
    printf "##fileformat=VCFv4.2\n"
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    printf "%s\t%s\t.\t%s\t%s\t.\t.\t.\n" \
      "$CHROM" "$pos" "$ref" "$alt"
  } > "$vcf"

  bgzip -f -c "$vcf" > "$vcfgz"
  tabix -f -p vcf "$vcfgz"

  local outfa="$OUTDIR/refs/ATCC19977_${label}.fasta"

  bcftools consensus \
    -s - \
    -f "$WT" \
    "$vcfgz" \
    > "$outfa"

  samtools faidx "$outfa"

  printf "%s\t%s\t%s\t%s\t%s\n" \
    "$label" \
    "$CHROM" \
    "$pos" \
    "$ref" \
    "$alt" \
    >> "$TRUTH"
}

# ------------------------------------------------------------
# rrl mutants
# ------------------------------------------------------------

for gp in "${RRL_SITES[@]}"; do
  make_mutant \
    "rrl_${gp}" \
    "$(gene_to_refpos "$RRL_START" "$gp")"
done

# ------------------------------------------------------------
# rrs mutants
# ------------------------------------------------------------

for gp in "${RRS_SITES[@]}"; do
  make_mutant \
    "rrs_${gp}" \
    "$(gene_to_refpos "$RRS_START" "$gp")"
done

# ------------------------------------------------------------
# erm41 SNP mutants
# ------------------------------------------------------------

# ATCC19977 erm41 position 19 is C.
# Simulate C19T.
pos19="$(gene_to_refpos "$ERM41_START" 19)"
make_mutant "erm41_C19T" "$pos19" "T"

# ATCC19977 erm41 position 28 is T.
# Simulate T28C.
pos28="$(gene_to_refpos "$ERM41_START" 28)"
make_mutant "erm41_T28C" "$pos28" "C"

# ------------------------------------------------------------
# Literature-described erm41 truncation
#
# Simulated genotype:
#
#   Δ64–65      = 2 bp deletion
#   Δ159–432    = 274 bp deletion
#
# Both coordinates are relative to the ORIGINAL erm41 gene in
# ATCC19977.
#
# Deletions are applied downstream -> upstream.
# ------------------------------------------------------------

make_erm41_truncation_fasta() {

  local label="$1"

  local small_start_ref
  local small_end_ref
  local large_start_ref
  local large_end_ref

  small_start_ref="$(
    gene_to_refpos "$ERM41_START" 64
  )"

  small_end_ref="$(
    gene_to_refpos "$ERM41_START" 65
  )"

  large_start_ref="$(
    gene_to_refpos "$ERM41_START" 159
  )"

  large_end_ref="$(
    gene_to_refpos "$ERM41_START" 432
  )"

  local outfa="$OUTDIR/refs/ATCC19977_${label}.fasta"

  echo
  echo "Creating literature-described erm41 truncation:"
  echo "  Δ64-65"
  echo "  Δ159-432"
  echo
  echo "Reference coordinates:"
  echo "  small: ${CHROM}:${small_start_ref}-${small_end_ref}"
  echo "  large: ${CHROM}:${large_start_ref}-${large_end_ref}"
  echo

  "$NTM_PYTHON" - \
    "$WT" \
    "$CHROM" \
    "$small_start_ref" \
    "$small_end_ref" \
    "$large_start_ref" \
    "$large_end_ref" \
    "$outfa" <<'PY'

import sys
from pathlib import Path

wt = Path(sys.argv[1])
chrom = sys.argv[2]

small_start = int(sys.argv[3])
small_end = int(sys.argv[4])

large_start = int(sys.argv[5])
large_end = int(sys.argv[6])

outfa = Path(sys.argv[7])

# Read FASTA
seqs = {}
order = []

name = None
buf = []

with wt.open() as f:
    for line in f:
        line = line.rstrip("\n")

        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)

            name = line[1:].split()[0]
            order.append(name)
            buf = []
        else:
            buf.append(line)

    if name is not None:
        seqs[name] = "".join(buf)

if chrom not in seqs:
    raise SystemExit(
        f"ERROR: contig {chrom} not found in {wt}"
    )

seq = seqs[chrom]

intervals = [
    ("small", small_start, small_end),
    ("large", large_start, large_end),
]

for label, start, end in intervals:
    if start < 1:
        raise SystemExit(
            f"ERROR: {label} deletion starts before sequence"
        )

    if end > len(seq):
        raise SystemExit(
            f"ERROR: {label} deletion extends beyond sequence "
            f"(end={end}, sequence length={len(seq)})"
        )

    if end < start:
        raise SystemExit(
            f"ERROR: invalid {label} deletion interval "
            f"{start}-{end}"
        )

# Delete downstream -> upstream.
deletions = [
    (small_start, small_end),
    (large_start, large_end),
]

deletions.sort(key=lambda x: x[0], reverse=True)

for start, end in deletions:
    start0 = start - 1
    end0 = end
    seq = seq[:start0] + seq[end0:]

seqs[chrom] = seq

# Write FASTA
with outfa.open("w") as out:
    for nm in order:
        out.write(f">{nm}\n")
        s = seqs[nm]

        for i in range(0, len(s), 60):
            out.write(s[i:i + 60] + "\n")

expected_removed = (
    (small_end - small_start + 1)
    +
    (large_end - large_start + 1)
)

print(
    f"Deleted {expected_removed} total bp "
    f"from {chrom}"
)

print(
    f"  small deletion: "
    f"{small_start}-{small_end} "
    f"({small_end-small_start+1} bp)"
)

print(
    f"  large deletion: "
    f"{large_start}-{large_end} "
    f"({large_end-large_start+1} bp)"
)

print(
    "Expected total deletion size: "
    f"{expected_removed} bp"
)
PY

  samtools faidx "$outfa"

  # ----------------------------------------------------------
  # SV truth table
  # ----------------------------------------------------------

  local svtruth="$OUTDIR/tmp/truth_sv.tsv"

  printf \
    "Dataset\tGene\tGene_DEL_start\tGene_DEL_end\tDEL_len\tCHROM\tRef_DEL_start\tRef_DEL_end\n" \
    > "$svtruth"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$label" \
    "erm41" \
    "64" \
    "65" \
    "2" \
    "$CHROM" \
    "$small_start_ref" \
    "$small_end_ref" \
    >> "$svtruth"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$label" \
    "erm41" \
    "159" \
    "432" \
    "274" \
    "$CHROM" \
    "$large_start_ref" \
    "$large_end_ref" \
    >> "$svtruth"

  echo
  echo "Created:"
  echo "  $outfa"
  echo
  echo "SV truth:"
  echo "  $svtruth"
}

make_erm41_truncation_fasta "erm41_truncated"

# ============================================================
# Additional erm-gene controls
# ============================================================

echo
echo "Creating additional erm-gene controls..."

# ------------------------------------------------------------
# ATCC35855 native erm39 control
#
# No sequence modification is performed. The genome itself
# contains erm39; this provides a native positive-control genome.
# ------------------------------------------------------------

ATCC35855_OUT="$OUTDIR/refs/ATCC35855_erm39.fasta"

cp -f "$ATCC35855_REF" "$ATCC35855_OUT"
samtools faidx "$ATCC35855_OUT"

echo "Created native erm39 control:"
echo "  $ATCC35855_OUT"

# ------------------------------------------------------------
# CCUG47445 synthetic erm55 insertion controls
#
# Complete chromosome:
#   CP007220.1
#
# Each of the three exact erm55 query sequences from
# references/nucleotide.fna is inserted AFTER nucleotide
# 2,500,000.
#
# This coordinate is intentionally synthetic. Its only purpose
# is to provide the same genomic background and insertion site
# for all three erm55 variants.
# ------------------------------------------------------------

CCUG_CONTIG="CP007220.1"
ERM55_INSERT_AFTER=2500000

make_erm55_insertion_fasta() {

  local erm_name="$1"
  local output_label="$2"

  local outfa="$OUTDIR/refs/CCUG47445_${output_label}.fasta"

  echo
  echo "Creating CCUG47445 erm55 control:"
  echo "  erm sequence:      $erm_name"
  echo "  chromosome:        $CCUG_CONTIG"
  echo "  insert after base: $ERM55_INSERT_AFTER"
  echo "  output:            $outfa"

  "$NTM_PYTHON" - \
    "$CCUG47445_REF" \
    "$NUCLEOTIDE_REF" \
    "$CCUG_CONTIG" \
    "$ERM55_INSERT_AFTER" \
    "$erm_name" \
    "$outfa" <<'PY'

import sys
from pathlib import Path

genome_path = Path(sys.argv[1])
query_path = Path(sys.argv[2])
target_contig = sys.argv[3]
insert_after = int(sys.argv[4])
erm_name = sys.argv[5]
out_path = Path(sys.argv[6])


def read_fasta(path):
    """
    Return:
      records = [
          {
              "id": first whitespace-delimited header token,
              "header": complete header text,
              "seq": sequence
          },
          ...
      ]
    """
    records = []

    header = None
    buf = []

    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    records.append(
                        {
                            "id": header.split()[0],
                            "header": header,
                            "seq": "".join(buf).upper(),
                        }
                    )

                header = line[1:]
                buf = []
            else:
                buf.append(line.strip())

        if header is not None:
            records.append(
                {
                    "id": header.split()[0],
                    "header": header,
                    "seq": "".join(buf).upper(),
                }
            )

    return records


genome_records = read_fasta(genome_path)
query_records = read_fasta(query_path)

query_lookup = {
    record["id"]: record["seq"]
    for record in query_records
}

if erm_name not in query_lookup:
    available = ", ".join(sorted(query_lookup))
    raise SystemExit(
        f"ERROR: {erm_name!r} not found in {query_path}. "
        f"Available IDs: {available}"
    )

erm_seq = query_lookup[erm_name]

if not erm_seq:
    raise SystemExit(
        f"ERROR: {erm_name} has an empty sequence"
    )

target_found = False

for record in genome_records:
    if record["id"] != target_contig:
        continue

    target_found = True
    seq = record["seq"]

    if insert_after < 1 or insert_after >= len(seq):
        raise SystemExit(
            f"ERROR: insertion coordinate {insert_after} is invalid "
            f"for {target_contig} length {len(seq)}"
        )

    # insert_after is a 1-based base coordinate.
    #
    # Example:
    #   insert_after = 2,500,000
    #
    # New sequence:
    #   original bases 1..2,500,000
    #   + erm55
    #   + original bases 2,500,001..end
    record["seq"] = (
        seq[:insert_after]
        + erm_seq
        + seq[insert_after:]
    )

    original_len = len(seq)
    new_len = len(record["seq"])

    expected_len = original_len + len(erm_seq)

    if new_len != expected_len:
        raise SystemExit(
            f"ERROR: insertion length check failed: "
            f"expected {expected_len}, observed {new_len}"
        )

    print(
        f"Inserted {erm_name}: "
        f"{len(erm_seq)} bp after "
        f"{target_contig}:{insert_after}"
    )

    print(
        f"Genome length: "
        f"{original_len} -> {new_len}"
    )

if not target_found:
    available = ", ".join(
        record["id"]
        for record in genome_records
    )

    raise SystemExit(
        f"ERROR: target contig {target_contig!r} not found in "
        f"{genome_path}. Available contigs: {available}"
    )

with out_path.open("w") as out:
    for record in genome_records:
        out.write(f">{record['header']}\n")

        seq = record["seq"]

        for i in range(0, len(seq), 60):
            out.write(seq[i:i + 60] + "\n")
PY

  samtools faidx "$outfa"

  echo "Created:"
  echo "  $outfa"
}

make_erm55_insertion_fasta \
  "erm55-plasmid" \
  "erm55_plasmid"

make_erm55_insertion_fasta \
  "erm55-transposon" \
  "erm55_transposon"

make_erm55_insertion_fasta \
  "erm55-chromosome" \
  "erm55_chromosome"

# ------------------------------------------------------------
# ATCC19977 Verification + Summary
#
# IMPORTANT:
# verify_sites.py is based on ATCC19977 coordinates.
# Therefore ONLY ATCC19977-derived FASTAs are included here.
#
# ATCC35855 and CCUG47445 controls are simulated for downstream
# erm gene detection and are intentionally excluded from this
# site-verification step.
# ------------------------------------------------------------

VERIFDIR="$OUTDIR/site_verification"

mkdir -p "$VERIFDIR"

VERIFY_SCRIPT="$REPO_ROOT/scripts/verify_sites.py"

if [[ ! -f "$VERIFY_SCRIPT" ]]; then
  echo "ERROR: Verification script not found: $VERIFY_SCRIPT" >&2
  exit 1
fi

# Clean old verification outputs so stale datasets from a
# previous simulation are not included.
rm -f \
  "$VERIFDIR"/*_site_verification.tsv \
  "$VERIFDIR"/site_verification_summary.tsv \
  2>/dev/null || true

echo
echo "Running ATCC19977 site verification..."

for fa in "$OUTDIR/refs/ATCC19977_"*.fasta; do

  [[ -f "$fa" ]] || continue

  "$NTM_PYTHON" \
    "$VERIFY_SCRIPT" \
    --fasta "$fa" \
    --outdir "$VERIFDIR" \
    >/dev/null

done

SUMMARY="$VERIFDIR/site_verification_summary.tsv"

"$NTM_PYTHON" - \
  "$VERIFDIR" \
  "$SUMMARY" <<'PY'

import re
import sys

from pathlib import Path

import pandas as pd


verifdir = Path(sys.argv[1]).resolve()
summary = Path(sys.argv[2]).resolve()


files = sorted(
    verifdir.glob("*_site_verification.tsv")
)

if not files:
    raise SystemExit(
        f"ERROR: No *_site_verification.tsv files found "
        f"in {verifdir}"
    )


def mutated_target_from_filename(
    tsv_path: Path
) -> str:

    # "<fasta minus .fasta>_site_verification.tsv"
    name = tsv_path.name.replace(
        "_site_verification.tsv",
        ""
    )

    # WT special case
    if (
        name.endswith("_WT")
        or name == "ATCC19977_WT"
    ):
        return "WT"

    # Strip leading ATCC19977_
    if name.startswith("ATCC19977_"):
        dataset = name[len("ATCC19977_"):]
    else:
        dataset = name

    # Keep truncation dataset name intact
    if dataset == "erm41_truncated":
        return "erm41_truncated"

    # erm41_C19T -> erm41_19
    # erm41_T28C -> erm41_28
    if dataset.startswith("erm41_"):

        nums = re.findall(
            r"(\d+)",
            dataset
        )

        if nums:
            return f"erm41_{nums[-1]}"

        return "erm41"

    # rrl_2269 and rrs_1458 already have desired names
    return dataset


out_rows = []

for f in files:

    df = pd.read_csv(
        f,
        sep="\t"
    )

    # Drop columns not needed in combined summary
    for col in [
        "fasta",
        "sequence_id",
        "ref_position_1based",
    ]:
        if col in df.columns:
            df = df.drop(
                columns=[col]
            )

    # Rename columns
    df = df.rename(
        columns={
            "gene_position_1based": "position",
            "expected_base": "atcc_strain",
        }
    )

    required = [
        "gene",
        "position",
        "atcc_strain",
        "observed_base",
        "match",
    ]

    missing = [
        c
        for c in required
        if c not in df.columns
    ]

    if missing:
        raise SystemExit(
            f"ERROR: {f.name} missing columns "
            f"{missing}. Found {list(df.columns)}"
        )

    df.insert(
        0,
        "mutated_target",
        mutated_target_from_filename(f)
    )

    df = df[
        [
            "mutated_target",
            "gene",
            "position",
            "atcc_strain",
            "observed_base",
            "match",
        ]
    ]

    out_rows.append(df)


out = pd.concat(
    out_rows,
    ignore_index=True
)

out.to_csv(
    summary,
    sep="\t",
    index=False
)

print(
    f"Wrote summary: {summary}"
)
PY

echo
echo "Done."
echo "Workdir:                 $WORKDIR"
echo "Input ATCC19977 FASTA:   $REF"
echo "Output dataset:          $OUTDIR"
echo "Simulation FASTAs:       $OUTDIR/refs/"
echo "SNP truth table:         $TRUTH"
echo "SV truth table:          $OUTDIR/tmp/truth_sv.tsv"
echo "Per-FASTA verification:  $VERIFDIR/*_site_verification.tsv"
echo "Summary verification:    $SUMMARY"
echo
echo "Additional controls:"
echo "  ATCC35855_erm39"
echo "  CCUG47445_erm55_plasmid"
echo "  CCUG47445_erm55_transposon"
echo "  CCUG47445_erm55_chromosome"
