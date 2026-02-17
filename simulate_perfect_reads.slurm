#!/usr/bin/env bash
#SBATCH --job-name=wgsim_perfect
#SBATCH --output=slurm_logs/wgsim_%A_%a.out
#SBATCH --error=slurm_logs/wgsim_%A_%a.err
#SBATCH --account=amr_services
#SBATCH --partition=standard
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -euo pipefail

# ------------------------------------------------------------
# simulate_perfect_reads.slurm
#
# Submits as a SLURM array job: one task per FASTA in
#   <workdir>/atcc_dataset/refs/*.fasta
#
# Output per FASTA:
#   <workdir>/atcc_dataset/reads/<dataset>/
#     <dataset>_R1.trim.fq.gz
#     <dataset>_R2.trim.fq.gz
#     <dataset>_SE.trim.fq.gz   (NOW: contains 1 FASTQ record copied from R1)
#
# dataset = FASTA basename minus ".fasta" and minus leading "ATCC19977_"
# ------------------------------------------------------------

: "${WORKDIR:?ERROR: Set WORKDIR env var, e.g. WORKDIR=/scratch/sgj4qr/mycobac_validation}"

PAIRS="${PAIRS:-1000000}"
READLEN="${READLEN:-150}"
INSMEAN="${INSMEAN:-300}"
INSSD="${INSSD:-30}"
SEED="${SEED:-12345}"
FORCE="${FORCE:-0}"

REFSDIR="${WORKDIR%/}/atcc_dataset/refs"
OUTROOT="${WORKDIR%/}/atcc_dataset/reads"
LOGDIR="${WORKDIR%/}/slurm_logs"

mkdir -p "$OUTROOT" "$LOGDIR"
mkdir -p "slurm_logs" 2>/dev/null || true  # in case you run from workdir

# Modules
if command -v module >/dev/null 2>&1; then
  module load samtools >/dev/null 2>&1 || true
fi

command -v wgsim >/dev/null 2>&1 || { echo "ERROR: wgsim not found (load samtools?)" >&2; exit 1; }
command -v gzip  >/dev/null 2>&1 || { echo "ERROR: gzip not found" >&2; exit 1; }

# Build FASTA list deterministically
mapfile -t FASTAS < <(ls -1 "$REFSDIR"/*.fasta | sort)
N="${#FASTAS[@]}"
[[ "$N" -gt 0 ]] || { echo "ERROR: No FASTAs found in $REFSDIR" >&2; exit 1; }

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
[[ "$TASK_ID" -ge 1 && "$TASK_ID" -le "$N" ]] || { echo "ERROR: SLURM_ARRAY_TASK_ID ($TASK_ID) out of range 1..$N" >&2; exit 1; }

fa="${FASTAS[$((TASK_ID-1))]}"

base="$(basename "$fa")"     # ATCC19977_erm41_C19T.fasta
name="${base%.fasta}"        # ATCC19977_erm41_C19T

if [[ "$name" == ATCC19977_* ]]; then
  dataset="${name#ATCC19977_}"   # erm41_C19T
else
  dataset="$name"
fi

outdir="$OUTROOT/$dataset"
mkdir -p "$outdir"

r1="$outdir/${dataset}_R1.trim.fq"
r2="$outdir/${dataset}_R2.trim.fq"
se="$outdir/${dataset}_SE.trim.fq"

r1gz="${r1}.gz"
r2gz="${r2}.gz"
segz="${se}.gz"

if [[ "$FORCE" -eq 0 ]] && [[ -s "$r1gz" && -s "$r2gz" && -s "$segz" ]]; then
  echo "Skipping (exists): $dataset"
  exit 0
fi

rm -f "$r1" "$r2" "$se" "$r1gz" "$r2gz" "$segz" 2>/dev/null || true

seed_this=$((SEED + TASK_ID))

echo "Task ${TASK_ID}/${N}"
echo "FASTA:    $fa"
echo "Dataset:  $dataset"
echo "Output:   $outdir"
echo "Params:   pairs=$PAIRS readlen=$READLEN insert=${INSMEAN}±${INSSD} seed=$seed_this force=$FORCE"
echo

wgsim \
  -N "$PAIRS" \
  -1 "$READLEN" -2 "$READLEN" \
  -d "$INSMEAN" -s "$INSSD" \
  -e 0 -r 0 -R 0 -X 0 \
  -S "$seed_this" \
  "$fa" "$r1" "$r2"

# Minimal change: make SE non-empty by copying the first FASTQ record from R1
awk 'NR<=4{print}' "$r1" > "$se"

gzip -f "$r1" "$r2" "$se"

[[ -s "$r1gz" && -s "$r2gz" && -s "$segz" ]] || { echo "ERROR: Missing gz outputs for $dataset" >&2; exit 1; }

echo "Done: $dataset"

