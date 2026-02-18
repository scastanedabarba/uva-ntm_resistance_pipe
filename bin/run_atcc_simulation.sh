#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# run_atcc_simulation.sh
#
# Wrapper to:
#   1) Run simulation/make_site_mutants.sh to generate WT + mutant FASTAs
#      (including deletion cases if added to make_site_mutants.sh)
#   2) Submit a SLURM array job to simulate PERFECT (error-free) reads
#      for each FASTA in <workdir>/atcc_dataset/refs/
#   3) Submit an additional small SLURM array job to simulate LOW-COVERAGE
#      reads for WT and one selected mutant (default: rrl_2270)
#   4) Submit a mixing job to create mixed-allele read datasets (5% and 50%)
#      by combining WT + mutant reads.
#
# Usage:
#   bash run_atcc_simulation.sh --workdir /scratch/sgj4qr/mycobac_validation
#
# Optional env overrides for read simulation:
#   PAIRS=1000000 READLEN=150 INSMEAN=300 INSSD=30 SEED=12345 FORCE=0
#
# Additional bundle overrides:
#   PAIRS_LOWCOV=100000   (default)  # low coverage pairs for WT + mutant
#   MIX_PAIRS=200000      (default)  # total pairs in mixed datasets
#   MUT_FOR_MIX=rrl_2270  (default)  # mutant dataset name (without ATCC19977_)
# ------------------------------------------------------------

usage() {
  cat <<'EOF'
run_atcc_simulation.sh

USAGE:
  bash run_atcc_simulation.sh --workdir <dir> [options]

REQUIRED:
  --workdir <dir>     Base working directory

OPTIONS:
  --chrom <name>      Contig name for make_site_mutants.sh (default: CU458896.1)
  --force_reads       Overwrite existing read files (sets FORCE=1 for SLURM jobs)
  -h, --help          Show help and exit

READ SIM PARAMS (via environment variables):
  PAIRS   (default 1000000)   # high coverage
  READLEN (default 150)
  INSMEAN (default 300)
  INSSD   (default 30)
  SEED    (default 12345)
  FORCE   (default 0)         # can also be set by --force_reads

MINIMAL BUNDLE PARAMS (via environment variables):
  PAIRS_LOWCOV (default 100000)   # low coverage for WT + MUT_FOR_MIX
  MIX_PAIRS    (default 200000)   # total pairs for mixed datasets
  MUT_FOR_MIX  (default rrl_2270) # dataset name in reads/ and refs/ (no ATCC19977_ prefix)
EOF
}

module load miniforge

WORKDIR=""
CHROM="CU458896.1"
FORCE_READS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir) WORKDIR="${2:-}"; shift 2;;
    --chrom) CHROM="${2:-}"; shift 2;;
    --force_reads) FORCE_READS=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "ERROR: Unknown option: $1" >&2; usage >&2; exit 1;;
  esac
done

[[ -z "$WORKDIR" ]] && { echo "ERROR: --workdir required" >&2; usage >&2; exit 1; }
WORKDIR="$(cd "$WORKDIR" && pwd)"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SIM_DIR="$REPO_ROOT/simulation"

SITE_MUTANTS="$SIM_DIR/make_site_mutants.sh"
SIM_SLURM="$SIM_DIR/simulate_reads_wgsim.slurm"
SIM_LIST_SLURM="$SIM_DIR/simulate_reads_from_list.slurm"
MIX_SLURM="$SIM_DIR/mix_reads.slurm"

[[ -f "$SITE_MUTANTS" ]] || { echo "ERROR: Missing: $SITE_MUTANTS" >&2; exit 1; }
[[ -f "$SIM_SLURM" ]] || { echo "ERROR: Missing SLURM script: $SIM_SLURM" >&2; exit 1; }

# Defaults for simulation (can be overridden via env)
PAIRS="${PAIRS:-1000000}"
READLEN="${READLEN:-150}"
INSMEAN="${INSMEAN:-300}"
INSSD="${INSSD:-30}"
SEED="${SEED:-12345}"
FORCE="${FORCE:-0}"

# Minimal bundle defaults (can be overridden via env)
PAIRS_LOWCOV="${PAIRS_LOWCOV:-100000}"
MIX_PAIRS="${MIX_PAIRS:-$PAIRS}"
MUT_FOR_MIX="${MUT_FOR_MIX:-rrl_2270}"

if [[ "$FORCE_READS" -eq 1 ]]; then
  FORCE=1
fi

echo "=== ATCC simulation pipeline ==="
echo "Workdir:      $WORKDIR"
echo "Chrom:        $CHROM"
echo
echo "HIGH-COV:  PAIRS=$PAIRS READLEN=$READLEN INSMEAN=$INSMEAN INSSD=$INSSD SEED=$SEED FORCE=$FORCE"
echo "LOW-COV:   PAIRS_LOWCOV=$PAIRS_LOWCOV (WT + ${MUT_FOR_MIX})"
echo "MIX:       MUT_FOR_MIX=$MUT_FOR_MIX MIX_PAIRS=$MIX_PAIRS (5% and 50%)"
echo

# ------------------------------------------------------------
# Step 1: generate mutant FASTAs + verification summary
# ------------------------------------------------------------
echo "== Step 1: Generating WT + mutant FASTAs =="
bash "$SITE_MUTANTS" --workdir "$WORKDIR" --chrom "$CHROM"

REFSDIR="$WORKDIR/atcc_dataset/refs"
[[ -d "$REFSDIR" ]] || { echo "ERROR: Expected refs dir not found: $REFSDIR" >&2; exit 1; }

# Count FASTAs
N_FASTA=$(ls -1 "$REFSDIR"/*.fasta 2>/dev/null | wc -l | awk '{print $1}')
[[ "$N_FASTA" -gt 0 ]] || { echo "ERROR: No FASTAs found in $REFSDIR" >&2; exit 1; }

echo
echo "Generated $N_FASTA FASTAs in: $REFSDIR"
echo

# ------------------------------------------------------------
# Step 2a: submit SLURM array job to simulate HIGH-COV reads
# ------------------------------------------------------------
mkdir -p "$WORKDIR/slurm_logs"

echo "== Step 2a: Submitting SLURM array job to simulate HIGH-COV reads =="
echo "Submitting array: 1-$N_FASTA"
echo "SLURM script: $SIM_SLURM"
echo

# Submit from WORKDIR so relative slurm_logs path works in the SBATCH header
cd "$WORKDIR"

JOB_SUBMIT_OUT=$(
  sbatch \
    --export=ALL,WORKDIR="$WORKDIR",PAIRS="$PAIRS",READLEN="$READLEN",INSMEAN="$INSMEAN",INSSD="$INSSD",SEED="$SEED",FORCE="$FORCE" \
    --array=1-"$N_FASTA" \
    "$SIM_SLURM"
)

echo "$JOB_SUBMIT_OUT"

# Extract job id (expects: "Submitted batch job <id>")
MAIN_JOBID="$(echo "$JOB_SUBMIT_OUT" | awk '{print $NF}')"
if ! [[ "$MAIN_JOBID" =~ ^[0-9]+$ ]]; then
  echo "WARNING: Could not parse main job id from: $JOB_SUBMIT_OUT" >&2
  MAIN_JOBID=""
fi

# ------------------------------------------------------------
# Step 2b: submit LOW-COV simulation for WT + one mutant
# ------------------------------------------------------------
if [[ -f "$SIM_LIST_SLURM" ]]; then
  echo
  echo "== Step 2b: Submitting LOW-COV simulation for WT + ${MUT_FOR_MIX} =="

  LOWCOV_LIST="$WORKDIR/atcc_dataset/tmp/lowcov_fastas.txt"
  mkdir -p "$(dirname "$LOWCOV_LIST")"

  WT_FA="$WORKDIR/atcc_dataset/refs/ATCC19977_WT.fasta"
  MUT_FA="$WORKDIR/atcc_dataset/refs/ATCC19977_${MUT_FOR_MIX}.fasta"

  [[ -f "$WT_FA" ]]  || { echo "ERROR: WT FASTA missing: $WT_FA" >&2; exit 1; }
  [[ -f "$MUT_FA" ]] || { echo "ERROR: Mutant FASTA missing: $MUT_FA (set MUT_FOR_MIX?)" >&2; exit 1; }

  {
    echo "$WT_FA"
    echo "$MUT_FA"
  } > "$LOWCOV_LIST"

  LOW_JOB_OUT=$(
    sbatch \
      --export=ALL,WORKDIR="$WORKDIR",FASTA_LIST="$LOWCOV_LIST",PAIRS="$PAIRS_LOWCOV",READLEN="$READLEN",INSMEAN="$INSMEAN",INSSD="$INSSD",SEED="$SEED",FORCE="$FORCE",SUFFIX="lowcov" \
      --array=1-2 \
      "$SIM_LIST_SLURM"
  )
  echo "$LOW_JOB_OUT"
else
  echo
  echo "NOTE: Skipping LOW-COV simulation (missing $SIM_LIST_SLURM)" >&2
fi

# ------------------------------------------------------------
# Step 3: submit MIXING job (depends on high-cov reads existing)
# ------------------------------------------------------------
if [[ -f "$MIX_SLURM" ]]; then
  echo
  echo "== Step 3: Submitting read-mixing job (5% and 50%) =="

  # If we have a job id, use dependency so we only mix after reads exist
  if [[ -n "$MAIN_JOBID" ]]; then
    MIX_OUT=$(
      sbatch \
        --dependency=afterok:"$MAIN_JOBID" \
        --export=ALL,WORKDIR="$WORKDIR",WT_DATASET="WT",MUT_DATASET="$MUT_FOR_MIX",MIX_PAIRS="$MIX_PAIRS",SEED="$SEED",FORCE="$FORCE" \
        "$MIX_SLURM"
    )
  else
    MIX_OUT=$(
      sbatch \
        --export=ALL,WORKDIR="$WORKDIR",WT_DATASET="WT",MUT_DATASET="$MUT_FOR_MIX",MIX_PAIRS="$MIX_PAIRS",SEED="$SEED",FORCE="$FORCE" \
        "$MIX_SLURM"
    )
  fi
  echo "$MIX_OUT"
else
  echo
  echo "NOTE: Skipping read-mixing (missing $MIX_SLURM)" >&2
fi

# ------------------------------------------------------------
# Minimal add: write 1-column linelist from refs/*.fasta (reads may not exist yet)
# ------------------------------------------------------------
LINELIST_OUT="$WORKDIR/atcc_dataset/linelist.tsv"

{
  echo -e "isolate"
  for fa in "$REFSDIR"/*.fasta; do
    base="$(basename "$fa")"
    name="${base%.fasta}"
    if [[ "$name" == ATCC19977_* ]]; then
      echo "${name#ATCC19977_}"
    else
      echo "$name"
    fi
  done | sort
} > "$LINELIST_OUT"

echo
echo "Wrote linelist: $LINELIST_OUT"
echo

echo "Done submitting jobs."
echo "Reads will be written to: $WORKDIR/atcc_dataset/reads/"
echo "Logs: $WORKDIR/slurm_logs/ (and/or slurm_logs/ relative to submission dir)"

