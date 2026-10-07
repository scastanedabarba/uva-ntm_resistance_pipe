#!/usr/bin/env bash
#SBATCH --job-name=myco_compile_interpret
#SBATCH --partition=standard
#SBATCH --account=amr_services_paid
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "USAGE:"
  echo "  sbatch 03_compile_and_interpret.slurm <outdir> <linelist> <targets_fasta> [ref_contig]"
  exit 2
fi

OUTDIR="$1"
LINELIST="$2"
TARGETS_FASTA="$3"
REF_CONTIG="${4:-CU458896.1}"

# Resolve python script from WORKDIR/scripts
: "${REPO_ROOT:?ERROR: REPO_ROOT env var not set (export from controller)}"
PY_SCRIPT="${REPO_ROOT%/}/scripts/compile_summary.py"
PY_SCRIPT2="${REPO_ROOT%/}/scripts/interpret_calls.py"

if [[ ! -s "$PY_SCRIPT" ]]; then
  echo "ERROR: Missing step3 python script: $PY_SCRIPT" >&2
  exit 2
fi

if [[ ! -s "$PY_SCRIPT2" ]]; then
  echo "ERROR: Missing step4 python script: $PY_SCRIPT2" >&2
  exit 2
fi

source "$REPO_ROOT/config/setup_environment.sh"
setup_python_env

PROVENANCE_HELPER="$REPO_ROOT/scripts/provenance.sh"
[[ -f "$PROVENANCE_HELPER" ]] || { echo "ERROR: Provenance helper not found: $PROVENANCE_HELPER" >&2; exit 1; }
source "$PROVENANCE_HELPER"

# Parse isolates from the same linelist used by the pipeline.
mapfile -t ISOLATES < <(awk -F'[,\t]' 'NR==1 && tolower($1) ~ /isolate/ {next} NF>=1 && $1!="" {gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1); print $1}' "$LINELIST")

for ISOLATE in "${ISOLATES[@]}"; do
  ISO_DIR="$OUTDIR/$ISOLATE"
  [[ -d "$ISO_DIR" ]] || continue
  provenance_init "$ISO_DIR"
  STEP3_LOG="$ISO_DIR/.provenance/step3.commands.log"
  : > "$STEP3_LOG"
  provenance_section "$STEP3_LOG" "STEP 3 - Compilation and Interpretation"
  provenance_log "$STEP3_LOG" "$NTM_PYTHON $PY_SCRIPT -o $OUTDIR -l $LINELIST -f $TARGETS_FASTA -r $REF_CONTIG"
  provenance_log "$STEP3_LOG" "$NTM_PYTHON $PY_SCRIPT2 --outdir $OUTDIR --linelist $LINELIST --ref-fasta ${REPO_ROOT%/}/references/ATCC19977.fasta"
done

"$NTM_PYTHON" "$PY_SCRIPT" \
  -o "$OUTDIR" \
  -l "$LINELIST" \
  -f "$TARGETS_FASTA" \
  -r "$REF_CONTIG"

"$NTM_PYTHON" "$PY_SCRIPT2" \
  --outdir "$OUTDIR" \
  --linelist "$LINELIST" \
  --ref-fasta "${REPO_ROOT%/}/references/ATCC19977.fasta"

# -----------------------------------------------------------------------------
# Finalize per-isolate provenance after compilation/interpretation succeeds.
# -----------------------------------------------------------------------------
PIPELINE_VERSION="$(git -C "$REPO_ROOT" describe --tags --always --dirty 2>/dev/null || echo UNAVAILABLE)"
GIT_COMMIT="$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo UNAVAILABLE)"
GIT_BRANCH="$(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD 2>/dev/null || echo UNAVAILABLE)"
if git -C "$REPO_ROOT" diff --quiet --ignore-submodules HEAD 2>/dev/null; then GIT_STATUS="clean"; else GIT_STATUS="modified"; fi

sha256_or_unavailable() {
  local f="$1"
  if [[ -s "$f" ]]; then sha256sum "$f" | awk '{print $1}'; else echo UNAVAILABLE; fi
}

for ISOLATE in "${ISOLATES[@]}"; do
  ISO_DIR="$OUTDIR/$ISOLATE"
  [[ -d "$ISO_DIR" ]] || continue
  PROV_DIR="$ISO_DIR/.provenance"

  {
    printf 'UVA NTM Resistance Pipeline - Command Log\n'
    printf 'Isolate: %s\n' "$ISOLATE"
    printf 'Generated: %s\n\n' "$(date '+%Y-%m-%d %H:%M:%S %z')"
    for f in "$PROV_DIR/step1.commands.log" "$PROV_DIR/step2.commands.log" "$PROV_DIR/step3.commands.log"; do
      if [[ -s "$f" ]]; then cat "$f"; else printf 'WARNING: Missing provenance component: %s\n\n' "$(basename "$f")"; fi
    done
  } > "$ISO_DIR/commands.log"

  {
    printf 'UVA NTM Resistance Pipeline\n'
    printf '===========================\n\n'
    printf 'PIPELINE\n'
    printf 'Version: %s\n' "$PIPELINE_VERSION"
    printf 'Git commit: %s\n' "$GIT_COMMIT"
    printf 'Git branch: %s\n' "$GIT_BRANCH"
    printf 'Git status: %s\n\n' "$GIT_STATUS"
    printf 'ANALYSIS\n'
    printf 'Isolate: %s\n' "$ISOLATE"
    printf 'Analysis completed: %s\n\n' "$(date '+%Y-%m-%d %H:%M:%S %z')"
    printf 'SOFTWARE\n'
    [[ -s "$PROV_DIR/step1.versions.txt" ]] && cat "$PROV_DIR/step1.versions.txt"
    [[ -s "$PROV_DIR/step2.versions.txt" ]] && cat "$PROV_DIR/step2.versions.txt"
    provenance_version_line /dev/stdout "Step 3 Python" "$NTM_PYTHON" --version
    printf '\nREFERENCES\n'
    printf 'Primary reference: %s\n' "${REPO_ROOT%/}/references/ATCC19977.fasta"
    printf 'Primary reference SHA256: %s\n' "$(sha256_or_unavailable "${REPO_ROOT%/}/references/ATCC19977.fasta")"
    printf 'Target FASTA: %s\n' "$TARGETS_FASTA"
    printf 'Target FASTA SHA256: %s\n' "$(sha256_or_unavailable "$TARGETS_FASTA")"
    printf 'erm39 reference: %s\n' "${REPO_ROOT%/}/references/ATCC35855.fasta"
    printf 'erm39 reference SHA256: %s\n' "$(sha256_or_unavailable "${REPO_ROOT%/}/references/ATCC35855.fasta")"
  } > "$ISO_DIR/pipeline_version.txt"
done

