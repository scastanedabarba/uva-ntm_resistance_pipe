#!/usr/bin/env bash
#SBATCH --job-name=myco_step3_compile
#SBATCH --partition=standard
#SBATCH --account=amr_services_paid
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "USAGE:"
  echo "  sbatch step3_compile.slurm <workdir> <outdir> <linelist> <targets_fasta> [only_isolate] [ref_contig]"
  exit 2
fi

WORKDIR="$1"
OUTDIR="$2"
LINELIST="$3"
TARGETS_FASTA="$4"
ONLY_ISO="${5:-}"
REF_CONTIG="${6:-CU458896.1}"

# Resolve python script from WORKDIR/scripts
PY_SCRIPT="${WORKDIR%/}/scripts/step3_compile.py"
PY_SCRIPT2="${WORKDIR%/}/scripts/step4_interpret.py"

if [[ ! -s "$PY_SCRIPT" ]]; then
  echo "ERROR: Missing step3 python script: $PY_SCRIPT" >&2
  exit 2
fi

if [[ ! -s "$PY_SCRIPT2" ]]; then
  echo "ERROR: Missing step4 python script: $PY_SCRIPT2" >&2
  exit 2
fi

module load miniforge
module load biopython

python3 "$PY_SCRIPT" \
  -o "$OUTDIR" \
  -l "$LINELIST" \
  -f "$TARGETS_FASTA" \
  ${ONLY_ISO:+-s "$ONLY_ISO"} \
  -r "$REF_CONTIG"

python3 "$PY_SCRIPT2" \
  --outdir "$OUTDIR" \
  --ref-fasta "${WORKDIR%/}/references/ATCC19977.fasta"
