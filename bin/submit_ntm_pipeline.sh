#!/usr/bin/env bash
set -euo pipefail

# Repo-relative paths
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKFLOW_DIR="$REPO_ROOT/workflow"
SCRIPTS_DIR="$REPO_ROOT/scripts"
REF_DIR="$REPO_ROOT/references"
SIM_DIR="$REPO_ROOT/simulation"

# ------------------------------------------------------------
# submit_ntm_pipeline.sh  (STEP1 + STEP2 + STEP3)
#
# Inputs:
#   1) LINELIST   (CSV or TSV)
#   2) OUTDIR     (pipeline output root; per-isolate dirs written here)
#
# Linelist formats:
#   - Production: 2 columns: isolate, run   (CSV or TSV)
#   - Simulated : 1 column : isolate only   (allowed ONLY with --simulated <DIR>)
#
# Read inputs:
#   Production (default):
#     /project/amr_services/qc/<run>/<isolate>/
#       <isolate>_R1.trim.fq.gz
#       <isolate>_R2.trim.fq.gz
#
#   Simulated (--simulated <DIR>):
#     <DIR>/atcc_dataset/reads/<isolate>/
#       <isolate>_R1.trim.fq.gz
#       <isolate>_R2.trim.fq.gz
#
# References (packaged in repo by default):
#   repo/references/ATCC19977.fasta
#   repo/references/nucleotide.fna
#
# Submits:
#   - Step 1 per isolate (assembly + BLAST)
#   - Step 2 per isolate (mapping + variants + consensus)
#   - Step 3 once (compile + interpret) after all Step1+Step2 jobs finish
# ------------------------------------------------------------

usage() {
  cat <<'EOF'
submit_ntm_pipeline.sh

USAGE:
  bash submit_ntm_pipeline.sh <linelist.(csv|tsv)> <outdir> [options]

REQUIRED:
  linelist   CSV/TSV with isolate in first column (run required unless --simulated)
  outdir     Output directory root for pipeline outputs

OPTIONS:
  --simulated <DIR>  Run in ATCC simulated mode using reads in:
                      <DIR>/atcc_dataset/reads/<isolate>/<isolate>_R{1,2}.trim.fq.gz
                      In this mode, a 1-column linelist is allowed (run will be set to "SIM").
  --isolate <ID>     Run only a single isolate (must be present in linelist)
  --partition <p>    SLURM partition (default: standard)
  --ref-fasta <p>    Override reference FASTA (default: repo/references/ATCC19977.fasta)
  -h, --help         Show this help and exit
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 2 ]]; then
  echo "ERROR: Missing required arguments." >&2
  usage >&2
  exit 1
fi

LINELIST="$1"
OUTDIR="$2"

# Defaults
ACCOUNT="amr_services_paid"
PARTITION="standard"
ONLY_ISOLATE=""
SIMULATED_DIR=""
REF_FASTA_OVERRIDE=""

shift 2
while [[ $# -gt 0 ]]; do
  case "$1" in
    --partition)
      PARTITION="${2:-}"
      [[ -z "$PARTITION" ]] && { echo "ERROR: --partition requires a value" >&2; exit 1; }
      shift 2
      ;;
    --isolate)
      ONLY_ISOLATE="${2:-}"
      [[ -z "$ONLY_ISOLATE" ]] && { echo "ERROR: --isolate requires a value" >&2; exit 1; }
      shift 2
      ;;
    --simulated)
      SIMULATED_DIR="${2:-}"
      [[ -z "$SIMULATED_DIR" ]] && { echo "ERROR: --simulated requires a value" >&2; exit 1; }
      shift 2
      ;;
    --ref-fasta)
      REF_FASTA_OVERRIDE="${2:-}"
      [[ -z "$REF_FASTA_OVERRIDE" ]] && { echo "ERROR: --ref-fasta requires a value" >&2; exit 1; }
      shift 2
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

[[ -f "$LINELIST" ]] || { echo "ERROR: Linelist not found: $LINELIST" >&2; exit 1; }

# Normalize paths
LINELIST="$(cd "$(dirname "$LINELIST")" && pwd)/$(basename "$LINELIST")"
mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

# Simulated mode: normalize and validate
READS_ROOT=""
if [[ -n "$SIMULATED_DIR" ]]; then
  SIMULATED_DIR="$(cd "$SIMULATED_DIR" && pwd)"
  [[ -d "$SIMULATED_DIR" ]] || { echo "ERROR: --simulated is not a directory: $SIMULATED_DIR" >&2; exit 1; }
  READS_ROOT="$SIMULATED_DIR/atcc_dataset/reads"
  [[ -d "$READS_ROOT" ]] || { echo "ERROR: Simulated reads directory not found: $READS_ROOT" >&2; exit 1; }
fi

mkdir -p "$OUTDIR"/{logs,status,summary}

STATUS_TSV="$OUTDIR/status/isolate_status.tsv"
echo -e "Isolate\tRun\tR1\tR2\tReadStatus\tStep1Job\tStep2Job" > "$STATUS_TSV"

# Workflow + reference inputs
STEP1_SCRIPT="$WORKFLOW_DIR/01_assembly_blast.slurm"
STEP2_SCRIPT="$WORKFLOW_DIR/02_map_call.slurm"
STEP3_SCRIPT="$WORKFLOW_DIR/03_compile_and_interpret.slurm"
TARGETS_FASTA="$REF_DIR/nucleotide.fna"

for f in "$STEP1_SCRIPT" "$STEP2_SCRIPT" "$STEP3_SCRIPT" "$TARGETS_FASTA"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing required file: $f" >&2; exit 1; }
done

# Reference FASTA (repo-packaged default)
if [[ -n "$REF_FASTA_OVERRIDE" ]]; then
  REF_FASTA="$(cd "$(dirname "$REF_FASTA_OVERRIDE")" && pwd)/$(basename "$REF_FASTA_OVERRIDE")"
else
  REF_FASTA="$REF_DIR/ATCC19977.fasta"
fi
[[ -s "$REF_FASTA" ]] || { echo "ERROR: Reference FASTA not found or empty: $REF_FASTA" >&2; exit 1; }

# ------------------------------------------------------------
# Index the reference ONCE here (avoids race conditions + bwa index segfaults)
# ------------------------------------------------------------
echo "[$(date)] Pre-indexing reference once (repo default unless --ref-fasta): $REF_FASTA" >&2

if command -v module >/dev/null 2>&1; then
  module use /project/amr_services/modulefiles/ >/dev/null 2>&1 || true
  module load bwa >/dev/null 2>&1 || true
  module load samtools >/dev/null 2>&1 || true
fi

command -v bwa >/dev/null 2>&1 || { echo "ERROR: bwa not found in PATH (module load bwa?)" >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH (module load samtools?)" >&2; exit 1; }

if [[ ! -s "${REF_FASTA}.bwt" ]]; then
  bwa index "$REF_FASTA" > "$OUTDIR/logs/ref_bwa_index.out" 2> "$OUTDIR/logs/ref_bwa_index.err"
fi
if [[ ! -s "${REF_FASTA}.fai" ]]; then
  samtools faidx "$REF_FASTA" > "$OUTDIR/logs/ref_faidx.out" 2> "$OUTDIR/logs/ref_faidx.err"
fi

# ------------------------------------------------------------
# Parse linelist:
# - Accept CSV or TSV
# - Header optional
# - Emits: isolate<TAB>run
#   * If run missing (1 column), emits isolate<TAB>"" (blank run)
# ------------------------------------------------------------
parse_linelist() {
  awk -v FS='[,\t]' '
    BEGIN{OFS="\t"}
    NR==1 {
      if (tolower($1) ~ /isolate/ || tolower($2) ~ /run/) next
    }
    NF>=1 && $1!="" {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
      run=""
      if (NF>=2) {
        run=$2
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", run)
      }
      print $1, run
    }
  ' "$1"
}

FOUND_ONLY_ISOLATE=0
SAW_ONECOL=0
SAW_TWOCOL=0
STEP3_DEPS=()

while IFS=$'\t' read -r ISOLATE RUN; do
  if [[ -n "$ONLY_ISOLATE" && "$ISOLATE" != "$ONLY_ISOLATE" ]]; then
    continue
  fi
  [[ -n "$ONLY_ISOLATE" ]] && FOUND_ONLY_ISOLATE=1

  if [[ -z "${RUN:-}" ]]; then
    SAW_ONECOL=1
  else
    SAW_TWOCOL=1
  fi

  # Enforce: 1-column linelist allowed ONLY with --simulated
  if [[ -z "${RUN:-}" && -z "$READS_ROOT" ]]; then
    echo "ERROR: Linelist row for isolate '$ISOLATE' has no RUN (1-column format), but --simulated was not provided." >&2
    echo "       Provide a 2-column linelist (isolate,run) OR re-run with --simulated <DIR>." >&2
    exit 1
  fi

  # In simulated mode, run defaults to SIM if missing
  if [[ -n "$READS_ROOT" && -z "${RUN:-}" ]]; then
    RUN="SIM"
  fi

  # Resolve read directory
  if [[ -n "$READS_ROOT" ]]; then
    READDIR="${READS_ROOT}/${ISOLATE}"
  else
    READDIR="/project/amr_services/qc/${RUN}/${ISOLATE}"
  fi

  R1="${READDIR}/${ISOLATE}_R1.trim.fq.gz"
  R2="${READDIR}/${ISOLATE}_R2.trim.fq.gz"

  ISO_OUT="$OUTDIR/$ISOLATE"
  mkdir -p "$ISO_OUT"/{assembly,blast,mapping,variants,logs}

  if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    echo "WARNING: Missing PE reads for isolate '$ISOLATE' (run='${RUN}')." >&2
    echo -e "${ISOLATE}\t${RUN}\t${R1}\t${R2}\tMISSING_READS\t\t" >> "$STATUS_TSV"
    continue
  fi

  STEP1_JOBID="$(
    sbatch --parsable \
      --job-name="ntm1_${ISOLATE}" \
      --account="$ACCOUNT" \
      --partition="$PARTITION" \
      --output="$OUTDIR/logs/step1_${ISOLATE}_%j.out" \
      --error="$OUTDIR/logs/step1_${ISOLATE}_%j.err" \
      --export=ALL,TARGET_FASTA="$TARGETS_FASTA" \
      "$STEP1_SCRIPT" "$ISOLATE" "$RUN" "$R1" "$R2" "$OUTDIR"
  )"

  STEP2_JOBID="$(
    sbatch --parsable \
      --job-name="ntm2_${ISOLATE}" \
      --account="$ACCOUNT" \
      --partition="$PARTITION" \
      --output="$OUTDIR/logs/step2_${ISOLATE}_%j.out" \
      --error="$OUTDIR/logs/step2_${ISOLATE}_%j.err" \
      --export=ALL,REF_FASTA="$REF_FASTA" \
      "$STEP2_SCRIPT" "$ISOLATE" "$RUN" "$R1" "$R2" "$OUTDIR"
  )"

  STEP3_DEPS+=("$STEP1_JOBID" "$STEP2_JOBID")
  echo -e "${ISOLATE}\t${RUN}\t${R1}\t${R2}\tOK\t${STEP1_JOBID}\t${STEP2_JOBID}" >> "$STATUS_TSV"

done < <(parse_linelist "$LINELIST")

if [[ -n "$ONLY_ISOLATE" && "$FOUND_ONLY_ISOLATE" -eq 0 ]]; then
  echo "ERROR: Isolate '$ONLY_ISOLATE' not found in linelist: $LINELIST" >&2
  exit 1
fi

# ------------------------------------------------------------
# Submit Step 3 once, after all Step1+Step2 jobs finish
# ------------------------------------------------------------
if [[ ${#STEP3_DEPS[@]} -gt 0 ]]; then
  DEP_STR="$(IFS=:; echo "${STEP3_DEPS[*]}")"
  STEP3_TARGETS="$REF_DIR/nucleotide.fna"

  [[ -s "$SCRIPTS_DIR/compile_summary.py" ]] || { echo "ERROR: Missing compile script: $SCRIPTS_DIR/compile_summary.py" >&2; exit 1; }
  [[ -s "$SCRIPTS_DIR/interpret_calls.py" ]] || { echo "ERROR: Missing interpret script: $SCRIPTS_DIR/interpret_calls.py" >&2; exit 1; }
  [[ -s "$STEP3_TARGETS" ]] || { echo "ERROR: Missing targets fasta: $STEP3_TARGETS" >&2; exit 1; }

  sbatch \
    --job-name="ntm3_compile" \
    --account="$ACCOUNT" \
    --partition="$PARTITION" \
    --dependency="afterany:${DEP_STR}" \
    --chdir="$OUTDIR" \
    --output="$OUTDIR/logs/step3_compile_%j.out" \
    --error="$OUTDIR/logs/step3_compile_%j.err" \
    --export=ALL,REPO_ROOT="$REPO_ROOT" \
    "$STEP3_SCRIPT" \
      "$OUTDIR" \
      "$LINELIST" \
      "$STEP3_TARGETS" \
      "CU458896.1"
else
  echo "No isolate jobs were submitted." >&2
fi

echo "Pipeline submission complete."
echo "Repo root:   $REPO_ROOT"
echo "Reference:   $REF_FASTA"
echo "Outdir:      $OUTDIR"
echo "Status table: $STATUS_TSV"
if [[ -n "$READS_ROOT" ]]; then
  echo "Simulated:   $SIMULATED_DIR"
  echo "Reads root:  $READS_ROOT"
fi
if [[ "$SAW_ONECOL" -eq 1 && -n "$READS_ROOT" ]]; then
  echo "Note: 1-column linelist detected; RUN defaulted to 'SIM' where missing."
fi
