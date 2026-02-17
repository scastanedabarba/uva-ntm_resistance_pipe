#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# site_mutants.sh  (Step 1 ONLY)
#
# Creates WT + per-site mutant FASTAs from ATCC19977 reference, then:
#   - runs site_verification.py on each generated FASTA (including WT)
#   - concatenates all verification outputs into one summary TSV with schema:
#       mutated_target, gene, position, atcc_strain, observed_base, match
#
# Inputs:
#   <workdir>/references/ATCC19977.fasta
#
# Outputs:
#   <workdir>/atcc_dataset/
#     refs/
#     vcfs/
#     tmp/
#     site_verification/
#       <fastaNameMinusDotFasta>_site_verification.tsv
#       site_verification_summary.tsv
# ------------------------------------------------------------

usage() {
  cat <<'EOF'
USAGE:
  bash site_mutants.sh --workdir <path> [options]

REQUIRED:
  --workdir <dir>         Working directory (contains references/)

OPTIONAL:
  --chrom <name>          Contig name (default: CU458896.1)
  --no_skip_same          Do NOT skip when ALT==REF
  -h, --help              Show help
EOF
}

WORKDIR=""
CHROM="CU458896.1"
NO_SKIP_SAME=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir) WORKDIR="${2:-}"; shift 2;;
    --chrom) CHROM="${2:-}"; shift 2;;
    --no_skip_same) NO_SKIP_SAME=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "ERROR: Unknown option $1" >&2; usage >&2; exit 1;;
  esac
done

[[ -z "$WORKDIR" ]] && { echo "ERROR: --workdir required" >&2; exit 1; }

REF="${WORKDIR%/}/references/ATCC19977.fasta"
OUTDIR="${WORKDIR%/}/atcc_dataset"

[[ -f "$REF" ]] || { echo "ERROR: Reference FASTA not found: $REF" >&2; exit 1; }

# -----------------------
# Load modules (best-effort)
# -----------------------
if command -v module >/dev/null 2>&1; then
  module load htslib   >/dev/null 2>&1 || true
  module load samtools >/dev/null 2>&1 || true
  module load bcftools >/dev/null 2>&1 || true
  module load python   >/dev/null 2>&1 || true
fi

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing $1" >&2; exit 1; }; }
need samtools
need bcftools
need bgzip
need tabix
need python3

mkdir -p "$OUTDIR"/{refs,vcfs,tmp}

# -----------------------
# Gene starts (1-based)
# -----------------------
RRL_START=1464208
RRS_START=1462398
ERM41_START=2345955

# -----------------------
# Sites of interest
# -----------------------
RRL_SITES=(2269 2270 2271 2281 2293)
RRS_SITES=(1373 1375 1376 1458)

gene_to_refpos() { echo $(( $1 + $2 - 1 )); }

# Index reference
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

# Confirm contig exists
if ! cut -f1 "${REF}.fai" | grep -Fxq "$CHROM"; then
  echo "ERROR: Contig '$CHROM' not found in FASTA index" >&2
  echo "Available contigs (first 20):" >&2
  cut -f1 "${REF}.fai" | head -n 20 >&2
  exit 1
fi

get_ref_base() {
  samtools faidx "$REF" "${CHROM}:$1-$1" | awk 'NR==2{print toupper($0)}'
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

# Copy WT
WT="$OUTDIR/refs/ATCC19977_WT.fasta"
cp -f "$REF" "$WT"
samtools faidx "$WT"

TRUTH="$OUTDIR/tmp/truth_table.tsv"
printf "Dataset\tCHROM\tPOS_ref\tREF\tALT\n" > "$TRUTH"

make_mutant() {
  local label="$1"
  local pos="$2"
  local forced_alt="${3:-}"

  local ref alt
  ref="$(get_ref_base "$pos")"
  [[ -z "$ref" ]] && { echo "WARNING: no REF at $pos"; return 0; }

  if [[ -n "$forced_alt" ]]; then
    alt="$forced_alt"
  else
    alt="$(flip_base "$ref")"
  fi

  if [[ "$alt" == "$ref" && "$NO_SKIP_SAME" -eq 0 ]]; then
    echo "WARNING: ALT==REF at ${CHROM}:${pos} for ${label}. Skipping." >&2
    return 0
  fi

  local vcf="$OUTDIR/vcfs/${label}.vcf"
  local vcfgz="${vcf}.gz"

  {
    printf "##fileformat=VCFv4.2\n"
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    printf "%s\t%s\t.\t%s\t%s\t.\t.\t.\n" "$CHROM" "$pos" "$ref" "$alt"
  } > "$vcf"

  bgzip -f -c "$vcf" > "$vcfgz"
  tabix -f -p vcf "$vcfgz"

  local outfa="$OUTDIR/refs/ATCC19977_${label}.fasta"
  bcftools consensus -s - -f "$WT" "$vcfgz" > "$outfa"
  samtools faidx "$outfa"

  printf "%s\t%s\t%s\t%s\t%s\n" "$label" "$CHROM" "$pos" "$ref" "$alt" >> "$TRUTH"
}

# rrl mutants
for gp in "${RRL_SITES[@]}"; do
  make_mutant "rrl_${gp}" "$(gene_to_refpos "$RRL_START" "$gp")"
done

# rrs mutants (includes 1458)
for gp in "${RRS_SITES[@]}"; do
  make_mutant "rrs_${gp}" "$(gene_to_refpos "$RRS_START" "$gp")"
done

# erm41 mutants: force opposite of ATCC alleles
pos19="$(gene_to_refpos "$ERM41_START" 19)"
make_mutant "erm41_C19T" "$pos19" "T"

pos28="$(gene_to_refpos "$ERM41_START" 28)"
make_mutant "erm41_T28C" "$pos28" "C"


# -----------------------
# erm41 deletion test cases (SV-ish)
# -----------------------
# Small deletion: 5 bp inside erm41 (gene pos 100..104)
# Large deletion: 300 bp inside erm41 (gene pos 50..349)
# These will shift downstream coordinates inside erm41, which is intentional for testing.

make_deletion_fasta() {
  local label="$1"
  local del_start_ref="$2"   # 1-based ref coordinate
  local del_len="$3"         # length to delete

  local outfa="$OUTDIR/refs/ATCC19977_${label}.fasta"

  python3 - "$WT" "$CHROM" "$del_start_ref" "$del_len" "$outfa" <<'PY'
import sys
from pathlib import Path

wt = Path(sys.argv[1])
chrom = sys.argv[2]
start = int(sys.argv[3])      # 1-based
length = int(sys.argv[4])
outfa = Path(sys.argv[5])

# Read FASTA (simple parser)
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
    raise SystemExit(f"ERROR: contig {chrom} not found in {wt}")

s = seqs[chrom]
i0 = start - 1
i1 = i0 + length
if i0 < 0 or i1 > len(s):
    raise SystemExit(f"ERROR: deletion {chrom}:{start}-{start+length-1} out of bounds (len={len(s)})")

seqs[chrom] = s[:i0] + s[i1:]

with outfa.open("w") as out:
    for nm in order:
        out.write(f">{nm}\n")
        seq = seqs[nm]
        for j in range(0, len(seq), 60):
            out.write(seq[j:j+60] + "\n")
PY

  samtools faidx "$outfa"

  # Record in a separate SV truth table (do NOT force into the SNP truth table schema)
  local svtruth="$OUTDIR/tmp/truth_sv.tsv"
  if [[ ! -f "$svtruth" ]]; then
    printf "Dataset\tCHROM\tDEL_start_ref\tDEL_len\n" > "$svtruth"
  fi
  printf "%s\t%s\t%s\t%s\n" "$label" "$CHROM" "$del_start_ref" "$del_len" >> "$svtruth"
}

# Compute deletion coordinates in reference space
del_small_start="$(gene_to_refpos "$ERM41_START" 100)"
make_deletion_fasta "erm41_del_small" "$del_small_start" 5

del_large_start="$(gene_to_refpos "$ERM41_START" 50)"
make_deletion_fasta "erm41_del_large" "$del_large_start" 300


# ------------------------------------------------------------
# Verification + Summary
# ------------------------------------------------------------
VERIFDIR="$OUTDIR/site_verification"
mkdir -p "$VERIFDIR"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERIFY_SCRIPT="$SCRIPT_DIR/site_verification.py"

if [[ ! -f "$VERIFY_SCRIPT" ]]; then
  echo "ERROR: Verification script not found: $VERIFY_SCRIPT" >&2
  exit 1
fi

# Clean old outputs to avoid stale/truncated filenames being included
rm -f "$VERIFDIR"/*_site_verification.tsv "$VERIFDIR"/site_verification_summary.tsv 2>/dev/null || true

echo
echo "Running site verification on generated FASTAs..."
for fa in "$OUTDIR/refs/"*.fasta; do
  python3 "$VERIFY_SCRIPT" --fasta "$fa" --outdir "$VERIFDIR" >/dev/null
done

SUMMARY="$VERIFDIR/site_verification_summary.tsv"

python3 - "$VERIFDIR" "$SUMMARY" <<'PY'
import re
import sys
from pathlib import Path
import pandas as pd

verifdir = Path(sys.argv[1]).resolve()
summary = Path(sys.argv[2]).resolve()

files = sorted(verifdir.glob("*_site_verification.tsv"))
if not files:
    raise SystemExit(f"ERROR: No *_site_verification.tsv files found in {verifdir}")

def mutated_target_from_filename(tsv_path: Path) -> str:
    # "<fasta minus .fasta>_site_verification.tsv"
    name = tsv_path.name.replace("_site_verification.tsv", "")

    # WT special case
    if name.endswith("_WT") or name == "ATCC19977_WT":
        return "WT"

    # Strip leading "ATCC19977_" if present
    dataset = name[len("ATCC19977_"):] if name.startswith("ATCC19977_") else name

    # erm41_C19T -> erm41_19 ; erm41_T28C -> erm41_28
    if dataset.startswith("erm41_"):
        nums = re.findall(r"(\d+)", dataset)
        if nums:
            return f"erm41_{nums[-1]}"  # LAST number is 19 or 28 (not 41)
        return "erm41"

    # rrl_2269, rrs_1458 already correct
    return dataset

out_rows = []
for f in files:
    df = pd.read_csv(f, sep="\t")

    # Drop columns
    for col in ["fasta", "sequence_id", "ref_position_1based"]:
        if col in df.columns:
            df = df.drop(columns=[col])

    # Rename columns
    df = df.rename(columns={
        "gene_position_1based": "position",
        "expected_base": "atcc_strain",
    })

    required = ["gene", "position", "atcc_strain", "observed_base", "match"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: {f.name} missing columns {missing}. Found {list(df.columns)}")

    df.insert(0, "mutated_target", mutated_target_from_filename(f))
    df = df[["mutated_target", "gene", "position", "atcc_strain", "observed_base", "match"]]
    out_rows.append(df)

out = pd.concat(out_rows, ignore_index=True)
out.to_csv(summary, sep="\t", index=False)
print(f"Wrote summary: {summary}")
PY

echo
echo "Done (Step 1)."
echo "Workdir:                 $WORKDIR"
echo "Input FASTA:             $REF"
echo "Output dataset:          $OUTDIR"
echo "Mutant FASTAs:           $OUTDIR/refs/"
echo "Truth table:             $TRUTH"
echo "Per-FASTA verification:  $VERIFDIR/*_site_verification.tsv"
echo "Summary verification:    $SUMMARY"
echo

