#!/usr/bin/env bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=amr_services_paid

set -euo pipefail

# ------------------------------------------------------------
# step2_map_call.slurm
#
# Step 2 (per isolate):
#   2a) Map PE reads to ATCC 19977 reference (bwa mem)
#   2b) Call variants across full target gene intervals (rrs, rrl, erm41)
#   2c) Emit per-gene consensus FASTA sequences for rrs/rrl/erm41
#
# NOTE:
#   - Reference indexing (bwa index + samtools faidx) is performed ONCE
#     by the controller. This script does NOT index the reference.
#   - This bcftools build does NOT support "-r" for consensus; we use
#     samtools faidx (region) piped into bcftools consensus instead.
# ------------------------------------------------------------

ISOLATE="${1:?Missing isolate}"
RUN="${2:?Missing run}"
R1FASTQ="${3:?Missing R1}"
R2FASTQ="${4:?Missing R2}"
PIPE_OUTDIR="${5:?Missing outdir}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# ---------------------------
# Modules
# ---------------------------
module use /project/amr_services/modulefiles/
module load bwa
module load samtools
module load bcftools

# ---------------------------
# Reference
# ---------------------------
# Controller should export REF_FASTA=<workdir>/references/ATCC19977.fasta
REF_FASTA="${REF_FASTA:-${PIPE_OUTDIR}/scripts/ATCC19977.fasta}"

# ---------------------------
# Target gene coordinates (ATCC 19977 reference numbering; 1-based inclusive)
# ---------------------------
RRS_START=1462398; RRS_END=1463901
RRL_START=1464208; RRL_END=1467319
ERM41_START=2345955; ERM41_END=2346476

# ---------------------------
# Outputs
# ---------------------------
ISO_DIR="${PIPE_OUTDIR}/${ISOLATE}"
MAP_DIR="${ISO_DIR}/mapping"
VAR_DIR="${ISO_DIR}/variants"
CONS_DIR="${VAR_DIR}/consensus"
LOG_DIR="${ISO_DIR}/logs"

mkdir -p "$MAP_DIR" "$VAR_DIR" "$CONS_DIR" "$LOG_DIR"

STATUS_TSV="${ISO_DIR}/status_step2.tsv"
echo -e "Isolate\tStep\tStatus\tMessage" > "$STATUS_TSV"

# ---------------------------
# Preflight
# ---------------------------
if [[ ! -s "$R1FASTQ" || ! -s "$R2FASTQ" ]]; then
  echo -e "${ISOLATE}\tSTEP2\tFAIL\tMissing R1/R2 input FASTQs" >> "$STATUS_TSV"
  exit 0
fi

if [[ ! -s "$REF_FASTA" ]]; then
  echo -e "${ISOLATE}\tSTEP2\tFAIL\tReference FASTA not found: $REF_FASTA" >> "$STATUS_TSV"
  exit 0
fi

# Require indexes to exist (controller should create them)
if [[ ! -s "${REF_FASTA}.bwt" ]]; then
  echo -e "${ISOLATE}\tBWA_INDEX\tFAIL\tMissing BWA index (.bwt). Controller should index reference first." >> "$STATUS_TSV"
  exit 0
fi
if [[ ! -s "${REF_FASTA}.fai" ]]; then
  echo -e "${ISOLATE}\tFAIDX\tFAIL\tMissing FASTA index (.fai). Controller should faidx reference first." >> "$STATUS_TSV"
  exit 0
fi

# Determine reference contig name from FASTA header (first token after '>')
REF_CONTIG="$(awk 'BEGIN{FS="[ \t]"} /^>/{gsub(/^>/,"",$1); print $1; exit}' "$REF_FASTA")"
if [[ -z "$REF_CONTIG" ]]; then
  echo -e "${ISOLATE}\tSTEP2\tFAIL\tCould not parse contig name from FASTA header" >> "$STATUS_TSV"
  exit 0
fi

# ---------------------------
# Build BED of target regions (0-based, half-open)
# ---------------------------
BED="${VAR_DIR}/targets.bed"
{
  printf "%s\t%d\t%d\trrs\n"   "$REF_CONTIG" $((RRS_START-1))   "$RRS_END"
  printf "%s\t%d\t%d\trrl\n"   "$REF_CONTIG" $((RRL_START-1))   "$RRL_END"
  printf "%s\t%d\t%d\term41\n" "$REF_CONTIG" $((ERM41_START-1)) "$ERM41_END"
} > "$BED"

# ---------------------------
# 2a) Mapping
# ---------------------------
echo "[$(date)] Mapping $ISOLATE with bwa mem" >&2

BAM="${MAP_DIR}/align.sorted.bam"

if bwa mem -t "$THREADS" "$REF_FASTA" "$R1FASTQ" "$R2FASTQ" 2> "${LOG_DIR}/bwa_mem.${ISOLATE}.err" \
  | samtools sort -@ "$THREADS" -o "$BAM" - 2> "${LOG_DIR}/samtools_sort.${ISOLATE}.err"; then
  :
else
  echo -e "${ISOLATE}\tMAP\tFAIL\tbwa mem / samtools sort failed (see logs)" >> "$STATUS_TSV"
  exit 0
fi

samtools index "$BAM" 2> "${LOG_DIR}/samtools_index.${ISOLATE}.err" || {
  echo -e "${ISOLATE}\tMAP\tFAIL\tsamtools index failed (see logs)" >> "$STATUS_TSV"
  exit 0
}

samtools flagstat "$BAM" > "${MAP_DIR}/flagstat.txt" 2> "${LOG_DIR}/flagstat.${ISOLATE}.err" || true
echo -e "${ISOLATE}\tMAP\tOK\tMapping complete" >> "$STATUS_TSV"

# ---------------------------
# 2b) Variant calling in targets (whole genes)
# ---------------------------
echo "[$(date)] Calling variants in target genes for $ISOLATE" >&2

VCF_GZ="${VAR_DIR}/targets.vcf.gz"

if bcftools mpileup \
    -f "$REF_FASTA" \
    -R "$BED" \
    -q 20 -Q 20 \
    -a DP,AD \
    -Ou \
    "$BAM" 2> "${LOG_DIR}/mpileup.${ISOLATE}.err" \
  | bcftools call -mv -Oz -o "$VCF_GZ" 2> "${LOG_DIR}/call.${ISOLATE}.err"; then
  :
else
  echo -e "${ISOLATE}\tVCF\tFAIL\tbcftools mpileup/call failed (see logs)" >> "$STATUS_TSV"
  exit 0
fi

bcftools index -t "$VCF_GZ" 2> "${LOG_DIR}/bcftools_index.${ISOLATE}.err" || {
  echo -e "${ISOLATE}\tVCF\tFAIL\tbcftools index failed (see logs)" >> "$STATUS_TSV"
  exit 0
}

DEPTH_TSV="${VAR_DIR}/targets_depth.tsv"
samtools depth -aa -b "$BED" "$BAM" > "$DEPTH_TSV" 2> "${LOG_DIR}/depth.${ISOLATE}.err" || true

RAW_VAR_TSV="${VAR_DIR}/targets_variants.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%AD]\n' "$VCF_GZ" \
  | awk -v iso="$ISOLATE" \
        -v rrs_s="$RRS_START" -v rrs_e="$RRS_END" \
        -v rrl_s="$RRL_START" -v rrl_e="$RRL_END" \
        -v erm_s="$ERM41_START" -v erm_e="$ERM41_END" \
        'BEGIN{FS=OFS="\t"; print "Isolate","Gene","CHROM","POS_ref","POS_gene","REF","ALT","QUAL","DP","AD"}
         {
           chrom=$1; pos=$2; ref=$3; alt=$4; qual=$5; dp=$6; ad=$7;
           gene="NA"; gpos="NA";
           if (pos>=rrs_s && pos<=rrs_e) { gene="rrs";  gpos=pos-rrs_s+1 }
           else if (pos>=rrl_s && pos<=rrl_e) { gene="rrl";  gpos=pos-rrl_s+1 }
           else if (pos>=erm_s && pos<=erm_e) { gene="erm41"; gpos=pos-erm_s+1 }
           print iso,gene,chrom,pos,gpos,ref,alt,qual,dp,ad
         }' > "$RAW_VAR_TSV"

echo -e "${ISOLATE}\tVCF\tOK\tVariant calling complete (targets only)" >> "$STATUS_TSV"

# ---------------------------
# 2c) Emit per-gene consensus FASTA (portable: samtools faidx | bcftools consensus)
# ---------------------------
echo "[$(date)] Writing per-gene consensus FASTA for $ISOLATE" >&2

write_consensus() {
  local gene="$1"
  local start="$2"
  local end="$3"
  local out_fa="$4"
  local region="${REF_CONTIG}:${start}-${end}"

  # Extract region FASTA, apply variants, rewrite header
  if samtools faidx "$REF_FASTA" "$region" 2> "${LOG_DIR}/faidx_region.${ISOLATE}.${gene}.err" \
    | bcftools consensus -f - "$VCF_GZ" 2> "${LOG_DIR}/consensus.${ISOLATE}.${gene}.err" \
    | awk -v h=">${ISOLATE}|${gene}|${region}" 'BEGIN{printed=0}
        /^>/{print h; printed=1; next}
        {print}
        END{if(!printed){print h}}' > "$out_fa"; then
    :
  else
    echo -e "${ISOLATE}\tCONSENSUS\tFAIL\tConsensus failed for ${gene} (see logs)" >> "$STATUS_TSV"
    return 1
  fi

  return 0
}

write_consensus "rrs"   "$RRS_START"   "$RRS_END"   "${CONS_DIR}/rrs.consensus.fasta"   || true
write_consensus "rrl"   "$RRL_START"   "$RRL_END"   "${CONS_DIR}/rrl.consensus.fasta"   || true
write_consensus "erm41" "$ERM41_START" "$ERM41_END" "${CONS_DIR}/erm41.consensus.fasta" || true

echo -e "${ISOLATE}\tCONSENSUS\tOK\tConsensus FASTAs written to ${CONS_DIR}" >> "$STATUS_TSV"

echo "[$(date)] Step 2 done for $ISOLATE" >&2

