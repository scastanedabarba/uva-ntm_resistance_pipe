#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# submit_ntm_pipeline.sh  (STEP1 + STEP2 + STEP3)
#
# Inputs:
#   1) LINELIST   (CSV or TSV)
#   2) WORKDIR    (must contain: references/ATCC19977.fasta)
#   3) OUTDIR     (pipeline output root; per-isolate dirs written here)
#
# Linelist formats:
#   - 2 columns: isolate, run   (CSV or TSV)   [always allowed]
#   - 1 column: isolate only    [allowed ONLY when --reads-root is provided]
#       (run will be set to "SIM")
#
# Read inputs:
#   Default (production):
#     /project/amr_services/qc/<run>/<isolate>/
#       <isolate>_R1.trim.fq.gz
#       <isolate>_R2.trim.fq.gz
#
#   With --reads-root:
#     <reads-root>/<isolate>/
#       <isolate>_R1.trim.fq.gz
#       <isolate>_R2.trim.fq.gz
#
# References (ALWAYS from WORKDIR):
#   WORKDIR/references/ATCC19977.fasta
#
# Submits:
#   - Step 1 per isolate (assembly + BLAST)
#   - Step 2 per isolate (mapping + variants + consensus)
#   - Step 3 once (compile) after all Step1+Step2 jobs finish
# ------------------------------------------------------------

usage() {
  cat <<'EOF'
submit_ntm_pipeline.sh

USAGE:
  bash submit_ntm_pipeline.sh <linelist.(csv|tsv)> <workdir> <outdir> [options]

REQUIRED:
  linelist   CSV/TSV with isolate in first column (run optional if --reads-root used)
  workdir    Base directory containing references/ATCC19977.fasta
  outdir     Output directory root for pipeline outputs

OPTIONS:
  --reads-root <DIR>  Override read root directory (useful for simulated reads)
  --isolate <ID>      Run only a single isolate (must be present in linelist)
  --partition <p>     SLURM partition (default: standard)
  -h, --help          Show this help and exit
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 3 ]]; then
  echo "ERROR: Missing required arguments." >&2
  usage >&2
  exit 1
fi

LINELIST="$1"
WORKDIR="$2"
OUTDIR="$3"

# Defaults
ACCOUNT="amr_services_paid"
PARTITION="standard"
ONLY_ISOLATE=""
READS_ROOT=""

shift 3
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
    --reads-root)
      READS_ROOT="${2:-}"
      [[ -z "$READS_ROOT" ]] && { echo "ERROR: --reads-root requires a value" >&2; exit 1; }
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
[[ -d "$WORKDIR" ]]  || { echo "ERROR: Workdir not found: $WORKDIR" >&2; exit 1; }

# Normalize WORKDIR and OUTDIR
WORKDIR="$(cd "$WORKDIR" && pwd)"
LINELIST="$(cd "$(dirname "$LINELIST")" && pwd)/$(basename "$LINELIST")"

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

# If reads-root provided, normalize to absolute path (and validate)
if [[ -n "$READS_ROOT" ]]; then
  READS_ROOT="$(cd "$READS_ROOT" && pwd)"
  [[ -d "$READS_ROOT" ]] || { echo "ERROR: --reads-root is not a directory: $READS_ROOT" >&2; exit 1; }
fi

mkdir -p "$OUTDIR"/{logs,status,summary,compiled}

STATUS_TSV="$OUTDIR/status/isolate_status.tsv"
echo -e "Isolate\tRun\tR1\tR2\tReadStatus\tStep1Job\tStep2Job" > "$STATUS_TSV"

# Resolve directory where this controller script lives
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

STEP1_SCRIPT="${SCRIPT_DIR}/step1_assembly_blast.slurm"
STEP2_SCRIPT="${SCRIPT_DIR}/step2_map_call.slurm"
STEP3_SCRIPT="${SCRIPT_DIR}/step3_compile.slurm"
TARGETS_FASTA="${SCRIPT_DIR}/nucleotide.fna"

for f in "$STEP1_SCRIPT" "$STEP2_SCRIPT" "$STEP3_SCRIPT" "$TARGETS_FASTA"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing required file: $f" >&2; exit 1; }
done

# Reference ALWAYS from workdir/references
REF_FASTA="${WORKDIR%/}/references/ATCC19977.fasta"
[[ -s "$REF_FASTA" ]] || { echo "ERROR: Reference FASTA not found: $REF_FASTA" >&2; exit 1; }

# ------------------------------------------------------------
# Index the reference ONCE here (avoids race conditions + bwa index segfaults)
# ------------------------------------------------------------
echo "[$(date)] Pre-indexing reference once: $REF_FASTA" >&2

if command -v module >/dev/null 2>&1; then
  module use /project/amr_services/modulefiles/ >/dev/null 2>&1 || true
  module load bwa >/dev/null 2>&1 || true
  module load samtools >/dev/null 2>&1 || true
fi

command -v bwa >/dev/null 2>&1 || { echo "ERROR: bwa not found in PATH (module load bwa?)" >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH (module load samtools?)" >&2; exit 1; }

# bwa index creates multiple files; use .bwt as sentinel
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

  # Enforce: 1-column linelist allowed ONLY with --reads-root
  if [[ -z "${RUN:-}" && -z "$READS_ROOT" ]]; then
    echo "ERROR: Linelist row for isolate '$ISOLATE' has no RUN (1-column format), but --reads-root was not provided." >&2
    echo "       Provide a 2-column linelist (isolate,run) OR re-run with --reads-root <DIR>." >&2
    exit 1
  fi

  # If reads-root is set and run is missing, use a dummy
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
  # keep consistent directory names (mapping not map)
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

  # Always use scripts from WORKDIR (avoids /var/spool/slurm path issues)
  STEP3_PY="${WORKDIR%/}/scripts/step3_compile.py"
  STEP3_TARGETS="${WORKDIR%/}/scripts/nucleotide.fna"

  [[ -s "$STEP3_PY" ]] || { echo "ERROR: Missing step3 python: $STEP3_PY" >&2; exit 1; }
  [[ -s "$STEP3_TARGETS" ]] || { echo "ERROR: Missing targets fasta: $STEP3_TARGETS" >&2; exit 1; }

  STEP3_ARGS=( "$WORKDIR" "$OUTDIR" "$LINELIST" "$STEP3_TARGETS" )
  if [[ -n "$ONLY_ISOLATE" ]]; then
    STEP3_ARGS+=( "$ONLY_ISOLATE" )
  else
    STEP3_ARGS+=( "" )
  fi
  STEP3_ARGS+=( "CU458896.1" )

  sbatch \
    --job-name="ntm3_compile" \
    --account="$ACCOUNT" \
    --partition="$PARTITION" \
    --dependency="afterany:${DEP_STR}" \
    --chdir="$WORKDIR" \
    --output="$OUTDIR/logs/step3_compile_%j.out" \
    --error="$OUTDIR/logs/step3_compile_%j.err" \
    "$STEP3_SCRIPT" "${STEP3_ARGS[@]}"

else
  echo "No isolate jobs were submitted." >&2
fi

echo "Pipeline submission complete."
echo "Workdir:     $WORKDIR"
echo "Reference:   $REF_FASTA"
echo "Outdir:      $OUTDIR"
echo "Status table: $STATUS_TSV"
if [[ -n "$READS_ROOT" ]]; then
  echo "Reads root:  $READS_ROOT"
fi
if [[ "$SAW_ONECOL" -eq 1 && -n "$READS_ROOT" ]]; then
  echo "Note: 1-column linelist detected; RUN defaulted to 'SIM' where missing."
fi

