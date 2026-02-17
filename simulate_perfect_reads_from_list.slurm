#!/usr/bin/env bash
#SBATCH --job-name=wgsim_perfect_list
#SBATCH --output=slurm_logs/wgsim_list_%A_%a.out
#SBATCH --error=slurm_logs/wgsim_list_%A_%a.err
#SBATCH --account=amr_services
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -euo pipefail

: "${WORKDIR:?ERROR: Set WORKDIR env var}"
: "${FASTA_LIST:?ERROR: Set FASTA_LIST env var (path to list of FASTA files)}"

PAIRS="${PAIRS:-100000}"
READLEN="${READLEN:-150}"
INSMEAN="${INSMEAN:-300}"
INSSD="${INSSD:-30}"
SEED="${SEED:-12345}"
FORCE="${FORCE:-0}"
SUFFIX="${SUFFIX:-lowcov}"   # appended to dataset name

OUTROOT="${WORKDIR%/}/atcc_dataset/reads"
LOGDIR="${WORKDIR%/}/slurm_logs"

mkdir -p "$OUTROOT" "$LOGDIR"
mkdir -p "slurm_logs" 2>/dev/null || true

if command -v module >/dev/null 2>&1; then
  module load samtools >/dev/null 2>&1 || true
fi

command -v wgsim >/dev/null 2>&1 || { echo "ERROR: wgsim not found" >&2; exit 1; }
command -v gzip  >/dev/null 2>&1 || { echo "ERROR: gzip not found" >&2; exit 1; }

mapfile -t FASTAS < "$FASTA_LIST"
N="${#FASTAS[@]}"
[[ "$N" -gt 0 ]] || { echo "ERROR: FASTA_LIST is empty: $FASTA_LIST" >&2; exit 1; }

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
[[ "$TASK_ID" -ge 1 && "$TASK_ID" -le "$N" ]] || { echo "ERROR: task id out of range" >&2; exit 1; }

fa="${FASTAS[$((TASK_ID-1))]}"
[[ -f "$fa" ]] || { echo "ERROR: FASTA not found: $fa" >&2; exit 1; }

base="$(basename "$fa")"
name="${base%.fasta}"
dataset="${name#ATCC19977_}"
dataset="${dataset}_${SUFFIX}"

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

seed_this=$((SEED + TASK_ID + 900000))

echo "FASTA: $fa"
echo "Dataset: $dataset"
echo "pairs=$PAIRS readlen=$READLEN insert=${INSMEAN}±${INSSD} seed=$seed_this"

wgsim \
  -N "$PAIRS" \
  -1 "$READLEN" -2 "$READLEN" \
  -d "$INSMEAN" -s "$INSSD" \
  -e 0 -r 0 -R 0 -X 0 \
  -S "$seed_this" \
  "$fa" "$r1" "$r2"

awk 'NR<=4{print}' "$r1" > "$se"
gzip -f "$r1" "$r2" "$se"

[[ -s "$r1gz" && -s "$r2gz" && -s "$segz" ]] || { echo "ERROR: missing outputs for $dataset" >&2; exit 1; }
echo "Done: $dataset"

