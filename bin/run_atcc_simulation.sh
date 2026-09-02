#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# run_atcc_simulation.sh
#
# Wrapper to:
#
#   1) Generate simulation FASTAs
#
#   2) Simulate HIGH-COVERAGE reads for every FASTA in:
#        <outdir>/atcc_dataset/refs/
#
#   3) Simulate LOW-COVERAGE reads for:
#        WT
#        rrl_2270
#
#   4) Create mixed read datasets:
#        rrl_2270_mix05
#        rrl_2270_mix50
#        rrs_1458_mix05
#        rrs_1458_mix50
#
#   5) After ALL read-generation jobs complete successfully,
#      create:
#
#        <outdir>/atcc_dataset/linelist.tsv
#
#      directly from valid read directories.
#
# This means read-derived controls such as low-coverage and
# mixed-allele datasets are automatically included in the
# downstream validation linelist.
#
# Usage:
#
#   bash run_atcc_simulation.sh \
#       --outdir /scratch/.../my_atcc_run
#
# Downstream:
#
#   bash bin/submit_ntm_pipeline.sh \
#       <linelist.tsv> \
#       <outdir> \
#       --simulated <outdir>
# ------------------------------------------------------------

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ------------------------------------------------------------
# Shared environment
# ------------------------------------------------------------

source "$REPO_ROOT/config/setup_environment.sh"

usage() {
  cat <<'EOF'
run_atcc_simulation.sh

USAGE:
  bash run_atcc_simulation.sh --outdir <dir> [options]

REQUIRED:
  --outdir <dir>
      Output directory where atcc_dataset/ will be created

OPTIONS:
  --chrom <name>
      ATCC19977 contig name for make_site_mutants.sh
      Default: CU458896.1

  --force_reads
      Overwrite existing read files
      Sets FORCE=1 for simulation jobs

  -h, --help
      Show help and exit
EOF
}

OUTDIR=""
CHROM="CU458896.1"
FORCE_READS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)
      OUTDIR="${2:-}"
      shift 2
      ;;
    --chrom)
      CHROM="${2:-}"
      shift 2
      ;;
    --force_reads)
      FORCE_READS=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

[[ -z "$OUTDIR" ]] && {
  echo "ERROR: --outdir required" >&2
  usage >&2
  exit 1
}

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

SIM_DIR="$REPO_ROOT/simulation"

SITE_MUTANTS="$SIM_DIR/make_site_mutants.sh"
SIM_SLURM="$SIM_DIR/simulate_reads_wgsim.slurm"
SIM_LIST_SLURM="$SIM_DIR/simulate_reads_from_list.slurm"
MIX_SLURM="$SIM_DIR/mix_reads.slurm"

[[ -f "$SITE_MUTANTS" ]] || {
  echo "ERROR: Missing: $SITE_MUTANTS" >&2
  exit 1
}

[[ -f "$SIM_SLURM" ]] || {
  echo "ERROR: Missing: $SIM_SLURM" >&2
  exit 1
}

[[ -f "$SIM_LIST_SLURM" ]] || {
  echo "ERROR: Missing: $SIM_LIST_SLURM" >&2
  exit 1
}

[[ -f "$MIX_SLURM" ]] || {
  echo "ERROR: Missing: $MIX_SLURM" >&2
  exit 1
}

# ------------------------------------------------------------
# Defaults
# ------------------------------------------------------------

PAIRS="${PAIRS:-1000000}"
READLEN="${READLEN:-150}"
INSMEAN="${INSMEAN:-300}"
INSSD="${INSSD:-30}"
SEED="${SEED:-12345}"
FORCE="${FORCE:-0}"

# Low coverage
PAIRS_LOWCOV="${PAIRS_LOWCOV:-100000}"
LOWCOV_MUTANT="${LOWCOV_MUTANT:-rrl_2270}"

# Mixed datasets
MIX_PAIRS="${MIX_PAIRS:-$PAIRS}"
RRL_MIX_MUTANT="${RRL_MIX_MUTANT:-rrl_2270}"
RRS_MIX_MUTANT="${RRS_MIX_MUTANT:-rrs_1458}"

if [[ "$FORCE_READS" -eq 1 ]]; then
  FORCE=1
fi

echo "=== NTM simulation pipeline ==="
echo
echo "Outdir: $OUTDIR"
echo "Chrom:  $CHROM"
echo
echo "HIGH-COV:"
echo "  PAIRS=$PAIRS"
echo "  READLEN=$READLEN"
echo "  INSMEAN=$INSMEAN"
echo "  INSSD=$INSSD"
echo "  SEED=$SEED"
echo "  FORCE=$FORCE"
echo
echo "LOW-COV:"
echo "  PAIRS_LOWCOV=$PAIRS_LOWCOV"
echo "  datasets:"
echo "    WT"
echo "    $LOWCOV_MUTANT"
echo
echo "MIXED:"
echo "  MIX_PAIRS=$MIX_PAIRS"
echo "  RRL mutant=$RRL_MIX_MUTANT"
echo "  RRS mutant=$RRS_MIX_MUTANT"
echo "  fractions=5% and 50%"
echo

# ------------------------------------------------------------
# Step 1
# Generate genome controls
# ------------------------------------------------------------

echo "== Step 1: Generating simulation FASTAs =="

bash "$SITE_MUTANTS" \
  --workdir "$OUTDIR" \
  --chrom "$CHROM"

REFSDIR="$OUTDIR/atcc_dataset/refs"

[[ -d "$REFSDIR" ]] || {
  echo "ERROR: Expected refs dir not found: $REFSDIR" >&2
  exit 1
}

mapfile -t REF_FASTAS < <(
  find "$REFSDIR" \
    -maxdepth 1 \
    -type f \
    -name '*.fasta' \
    | sort
)

N_FASTA="${#REF_FASTAS[@]}"

[[ "$N_FASTA" -gt 0 ]] || {
  echo "ERROR: No FASTAs found in $REFSDIR" >&2
  exit 1
}

echo
echo "Generated $N_FASTA simulation FASTAs:"
printf '  %s\n' "${REF_FASTAS[@]}"
echo

# ------------------------------------------------------------
# SLURM setup
# ------------------------------------------------------------

mkdir -p "$OUTDIR/slurm_logs"

cd "$OUTDIR"

# Track all terminal read-generation job IDs.
#
# The final linelist-generation job will depend on these jobs.
FINAL_DEPENDENCIES=()

# ------------------------------------------------------------
# Step 2a
# High-coverage simulation
# ------------------------------------------------------------

echo "== Step 2a: Submitting HIGH-COVERAGE read simulation =="

echo "Array: 1-$N_FASTA"
echo "SLURM script: $SIM_SLURM"
echo

JOB_SUBMIT_OUT=$(
  sbatch \
    --export=ALL,REPO_ROOT="$REPO_ROOT",WORKDIR="$OUTDIR",PAIRS="$PAIRS",READLEN="$READLEN",INSMEAN="$INSMEAN",INSSD="$INSSD",SEED="$SEED",FORCE="$FORCE" \
    --array=1-"$N_FASTA" \
    "$SIM_SLURM"
)

echo "$JOB_SUBMIT_OUT"

MAIN_JOBID="$(
  echo "$JOB_SUBMIT_OUT" |
  awk '{print $NF}'
)"

if ! [[ "$MAIN_JOBID" =~ ^[0-9]+$ ]]; then
  echo "ERROR: Could not parse high-coverage job ID from:" >&2
  echo "$JOB_SUBMIT_OUT" >&2
  exit 1
fi

FINAL_DEPENDENCIES+=("$MAIN_JOBID")

echo "High-coverage job ID: $MAIN_JOBID"

# ------------------------------------------------------------
# Step 2b
# Low-coverage simulation
#
# Keep:
#   WT_lowcov
#   rrl_2270_lowcov
#
# or another mutant if LOWCOV_MUTANT is overridden.
# ------------------------------------------------------------

echo
echo "== Step 2b: Submitting LOW-COVERAGE simulations =="

LOWCOV_LIST="$OUTDIR/atcc_dataset/tmp/lowcov_fastas.txt"

mkdir -p "$(dirname "$LOWCOV_LIST")"

WT_FA="$OUTDIR/atcc_dataset/refs/ATCC19977_WT.fasta"
LOWCOV_MUT_FA="$OUTDIR/atcc_dataset/refs/ATCC19977_${LOWCOV_MUTANT}.fasta"

[[ -f "$WT_FA" ]] || {
  echo "ERROR: WT FASTA missing: $WT_FA" >&2
  exit 1
}

[[ -f "$LOWCOV_MUT_FA" ]] || {
  echo "ERROR: Low-coverage mutant FASTA missing:" >&2
  echo "  $LOWCOV_MUT_FA" >&2
  exit 1
}

{
  echo "$WT_FA"
  echo "$LOWCOV_MUT_FA"
} > "$LOWCOV_LIST"

LOW_JOB_OUT=$(
  sbatch \
    --export=ALL,REPO_ROOT="$REPO_ROOT",WORKDIR="$OUTDIR",FASTA_LIST="$LOWCOV_LIST",PAIRS="$PAIRS_LOWCOV",READLEN="$READLEN",INSMEAN="$INSMEAN",INSSD="$INSSD",SEED="$SEED",FORCE="$FORCE",SUFFIX="lowcov" \
    --array=1-2 \
    "$SIM_LIST_SLURM"
)

echo "$LOW_JOB_OUT"

LOW_JOBID="$(
  echo "$LOW_JOB_OUT" |
  awk '{print $NF}'
)"

if ! [[ "$LOW_JOBID" =~ ^[0-9]+$ ]]; then
  echo "ERROR: Could not parse low-coverage job ID from:" >&2
  echo "$LOW_JOB_OUT" >&2
  exit 1
fi

FINAL_DEPENDENCIES+=("$LOW_JOBID")

echo "Low-coverage job ID: $LOW_JOBID"

# ------------------------------------------------------------
# Step 3a
# rrl mixed-allele simulations
#
# Depends on high-coverage WT + mutant reads.
# ------------------------------------------------------------

echo
echo "== Step 3a: Submitting rrl mixed-read simulations =="

RRL_MIX_OUT=$(
  sbatch \
    --dependency=afterok:"$MAIN_JOBID" \
    --export=ALL,REPO_ROOT="$REPO_ROOT",WORKDIR="$OUTDIR",WT_DATASET="WT",MUT_DATASET="$RRL_MIX_MUTANT",MIX_PAIRS="$MIX_PAIRS",SEED="$SEED",FORCE="$FORCE" \
    "$MIX_SLURM"
)

echo "$RRL_MIX_OUT"

RRL_MIX_JOBID="$(
  echo "$RRL_MIX_OUT" |
  awk '{print $NF}'
)"

if ! [[ "$RRL_MIX_JOBID" =~ ^[0-9]+$ ]]; then
  echo "ERROR: Could not parse rrl mixing job ID from:" >&2
  echo "$RRL_MIX_OUT" >&2
  exit 1
fi

FINAL_DEPENDENCIES+=("$RRL_MIX_JOBID")

echo "rrl mixing job ID: $RRL_MIX_JOBID"

# ------------------------------------------------------------
# Step 3b
# rrs mixed-allele simulations
#
# Generates:
#   rrs_1458_mix05
#   rrs_1458_mix50
# ------------------------------------------------------------

echo
echo "== Step 3b: Submitting rrs mixed-read simulations =="

RRS_MIX_OUT=$(
  sbatch \
    --dependency=afterok:"$MAIN_JOBID" \
    --export=ALL,REPO_ROOT="$REPO_ROOT",WORKDIR="$OUTDIR",WT_DATASET="WT",MUT_DATASET="$RRS_MIX_MUTANT",MIX_PAIRS="$MIX_PAIRS",SEED="$SEED",FORCE="$FORCE" \
    "$MIX_SLURM"
)

echo "$RRS_MIX_OUT"

RRS_MIX_JOBID="$(
  echo "$RRS_MIX_OUT" |
  awk '{print $NF}'
)"

if ! [[ "$RRS_MIX_JOBID" =~ ^[0-9]+$ ]]; then
  echo "ERROR: Could not parse rrs mixing job ID from:" >&2
  echo "$RRS_MIX_OUT" >&2
  exit 1
fi

FINAL_DEPENDENCIES+=("$RRS_MIX_JOBID")

echo "rrs mixing job ID: $RRS_MIX_JOBID"

# ------------------------------------------------------------
# Step 4
# Finalize linelist AFTER all read-generation jobs finish.
#
# The reads directory is the source of truth.
#
# A dataset is included only if both:
#
#   <dataset>_R1.trim.fq.gz
#   <dataset>_R2.trim.fq.gz
#
# exist and are non-empty.
# ------------------------------------------------------------

echo
echo "== Step 4: Submitting final linelist-generation job =="

DEPENDENCY_STRING="$(
  IFS=:
  echo "${FINAL_DEPENDENCIES[*]}"
)"

FINALIZE_CMD=$(cat <<'EOF'
set -euo pipefail

READSROOT="${WORKDIR%/}/atcc_dataset/reads"
LINELIST="${WORKDIR%/}/atcc_dataset/linelist.tsv"

[[ -d "$READSROOT" ]] || {
    echo "ERROR: reads directory not found: $READSROOT" >&2
    exit 1
}

TMP="${LINELIST}.tmp"

{
    echo "isolate"

    for dir in "$READSROOT"/*; do

        [[ -d "$dir" ]] || continue

        dataset="$(basename "$dir")"

        r1="$dir/${dataset}_R1.trim.fq.gz"
        r2="$dir/${dataset}_R2.trim.fq.gz"

        if [[ -s "$r1" && -s "$r2" ]]; then
            echo "$dataset"
        else
            echo "WARNING: excluding incomplete dataset: $dataset" >&2
        fi

    done | sort

} > "$TMP"

mv "$TMP" "$LINELIST"

echo
echo "Wrote linelist:"
echo "  $LINELIST"
echo

echo "Datasets:"
tail -n +2 "$LINELIST" | nl -ba

echo
echo -n "Total datasets: "
tail -n +2 "$LINELIST" | wc -l
EOF
)

LINELIST_JOB_OUT=$(
  sbatch \
    --job-name=sim_linelist \
    --output="$OUTDIR/slurm_logs/sim_linelist_%j.out" \
    --error="$OUTDIR/slurm_logs/sim_linelist_%j.err" \
    --account=amr_services_paid \
    --partition=standard \
    --time=00:10:00 \
    --cpus-per-task=1 \
    --mem=512M \
    --dependency=afterok:"$DEPENDENCY_STRING" \
    --export=ALL,WORKDIR="$OUTDIR" \
    --wrap="$FINALIZE_CMD"
)

echo "$LINELIST_JOB_OUT"

LINELIST_JOBID="$(
  echo "$LINELIST_JOB_OUT" |
  awk '{print $NF}'
)"

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Simulation jobs submitted"
echo "============================================================"
echo
echo "High coverage:"
echo "  $MAIN_JOBID"
echo
echo "Low coverage:"
echo "  $LOW_JOBID"
echo
echo "rrl mixtures:"
echo "  $RRL_MIX_JOBID"
echo
echo "rrs mixtures:"
echo "  $RRS_MIX_JOBID"
echo
echo "Final linelist:"
echo "  $LINELIST_JOBID"
echo
echo "Reads:"
echo "  $OUTDIR/atcc_dataset/reads/"
echo
echo "Linelist will be written after successful read simulation:"
echo "  $OUTDIR/atcc_dataset/linelist.tsv"
echo
echo "Logs:"
echo "  $OUTDIR/slurm_logs/"
echo
echo "Monitor with:"
echo "  squeue -u \$USER"
echo
echo "After the linelist job finishes, run:"
echo
echo "  bash $REPO_ROOT/bin/submit_ntm_pipeline.sh \\"
echo "      $OUTDIR/atcc_dataset/linelist.tsv \\"
echo "      $OUTDIR \\"
echo "      --simulated $OUTDIR"
echo
