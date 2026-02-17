#!/usr/bin/env bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail

# ------------------------------------------------------------
# step1_assembly_blast.slurm
#
# Step 1 (per isolate):
#   1) SPAdes assembly (PE + optional -s singletons)
#   2) BLAST: target genes (query) vs assembled contigs (db)
#      - Writes blast_raw.tsv (all hits)
#      - Writes blast_filtered.tsv (pident>=90 AND qcovs>=90)
#
# No interpretation is performed here.
# Interpretation will be done in Step 3 by parsing the outputs.
#
# Inputs (from controller):
#   $1 ISOLATE
#   $2 RUN
#   $3 R1FASTQ
#   $4 R2FASTQ
#   $5 OUTDIR (pipeline output root)
#
# Optional env var:
#   TARGET_FASTA=/path/to/scripts/nucleotide.fna
#     (defaults to <OUTDIR>/scripts/nucleotide.fna)
# ------------------------------------------------------------

ISOLATE="${1:?Missing isolate}"
RUN="${2:?Missing run}"
R1FASTQ="${3:?Missing R1}"
R2FASTQ="${4:?Missing R2}"
PIPE_OUTDIR="${5:?Missing outdir}"

# ---------------------------
# Modules
# ---------------------------
module use /project/amr_services/modulefiles/
module load spades/4.2.0
module load blast/2.17.0

# ---------------------------
# Outputs
# ---------------------------
ISO_DIR="${PIPE_OUTDIR}/${ISOLATE}"
ASM_DIR="${ISO_DIR}/assembly"
BLAST_DIR="${ISO_DIR}/blast"
LOG_DIR="${ISO_DIR}/logs"

mkdir -p "$ASM_DIR" "$BLAST_DIR" "$LOG_DIR"

STATUS_TSV="${ISO_DIR}/status_step1.tsv"
echo -e "Isolate\tStep\tStatus\tMessage" > "$STATUS_TSV"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Keep -s as requested; use the canonical location for SE reads
READDIR="$(dirname "$R1FASTQ")"
SEFASTQ="${READDIR}/${ISOLATE}_SE.trim.fq.gz"

# Use your FASTA exactly as-is
TARGET_FASTA="${TARGET_FASTA:-${PIPE_OUTDIR}/scripts/nucleotide.fna}"

# BLAST thresholds (do NOT interpret here; only filter)
MIN_PID=90
MIN_QCOV=90

# ---------------------------
# Preflight
# ---------------------------
if [[ ! -s "$R1FASTQ" || ! -s "$R2FASTQ" ]]; then
  echo -e "${ISOLATE}\tSTEP1\tFAIL\tMissing R1/R2 input FASTQs" >> "$STATUS_TSV"
  exit 0
fi

if [[ ! -s "$SEFASTQ" ]]; then
  echo "NOTE: SE file not found (will omit -s): $SEFASTQ" >&2
fi

if [[ ! -s "$TARGET_FASTA" ]]; then
  echo -e "${ISOLATE}\tSTEP1\tFAIL\tTarget FASTA not found: $TARGET_FASTA" >> "$STATUS_TSV"
  exit 0
fi

# ---------------------------
# 1) SPAdes assembly
# ---------------------------
echo "[$(date)] Running SPAdes for $ISOLATE"

SPADES_CMD=(spades.py --isolate --cov-cutoff auto -1 "$R1FASTQ" -2 "$R2FASTQ")
if [[ -s "$SEFASTQ" ]]; then
  SPADES_CMD+=(-s "$SEFASTQ")
fi
SPADES_CMD+=(-t "$THREADS" -o "$ASM_DIR")

if "${SPADES_CMD[@]}" > "${LOG_DIR}/spades.${ISOLATE}.out" 2> "${LOG_DIR}/spades.${ISOLATE}.err"; then
  echo -e "${ISOLATE}\tSPADES\tOK\tAssembly completed" >> "$STATUS_TSV"
else
  echo -e "${ISOLATE}\tSPADES\tFAIL\tSPAdes failed (see logs)" >> "$STATUS_TSV"
  exit 0
fi

CONTIGS="${ASM_DIR}/contigs.fasta"
if [[ ! -s "$CONTIGS" ]]; then
  echo -e "${ISOLATE}\tSPADES\tFAIL\tcontigs.fasta missing/empty" >> "$STATUS_TSV"
  exit 0
fi

# ---------------------------
# 2) BLAST (targets as query vs contigs as db)
# ---------------------------
echo "[$(date)] Running BLAST for $ISOLATE"

DB_PREFIX="${BLAST_DIR}/contigs_db"

makeblastdb -in "$CONTIGS" -dbtype nucl -out "$DB_PREFIX" \
  > "${LOG_DIR}/makeblastdb.${ISOLATE}.out" 2> "${LOG_DIR}/makeblastdb.${ISOLATE}.err" || {
    echo -e "${ISOLATE}\tBLASTDB\tFAIL\tmakeblastdb failed (see logs)" >> "$STATUS_TSV"
    exit 0
  }

RAW_TSV="${BLAST_DIR}/blast_raw.tsv"
FILTERED_TSV="${BLAST_DIR}/blast_filtered.tsv"

# Columns:
# qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs
blastn -query "$TARGET_FASTA" -db "$DB_PREFIX" \
  -task blastn \
  -max_target_seqs 50 \
  -evalue 1e-10 \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs" \
  > "$RAW_TSV" 2> "${LOG_DIR}/blastn.${ISOLATE}.err" || {
    echo -e "${ISOLATE}\tBLAST\tFAIL\tblastn failed (see logs)" >> "$STATUS_TSV"
    exit 0
  }

# Filter: pident >= 90 AND qcovs >= 90
awk -v pid="$MIN_PID" -v qcov="$MIN_QCOV" 'BEGIN{FS=OFS="\t"} ($3+0)>=pid && ($13+0)>=qcov {print}' \
  "$RAW_TSV" > "$FILTERED_TSV"

echo -e "${ISOLATE}\tBLAST\tOK\tBLAST completed; filtered with pident>=${MIN_PID} and qcovs>=${MIN_QCOV}" >> "$STATUS_TSV"

echo "[$(date)] Step 1 done for $ISOLATE"

