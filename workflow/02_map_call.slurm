#!/usr/bin/env bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --partition=standard
#SBATCH --account=amr_services_paid

set -euo pipefail

# ------------------------------------------------------------
# 02_map_call.slurm
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
# Environment
# ---------------------------
: "${REPO_ROOT:?ERROR: REPO_ROOT env var not set (export from controller)}"

ENV_CONFIG="$REPO_ROOT/config/setup_environment.sh"

[[ -f "$ENV_CONFIG" ]] || {
    echo "ERROR: Environment config not found: $ENV_CONFIG" >&2
    exit 1
}
source "$ENV_CONFIG"
setup_mapping_env
setup_python_env

PROVENANCE_HELPER="$REPO_ROOT/scripts/provenance.sh"
[[ -f "$PROVENANCE_HELPER" ]] || { echo "ERROR: Provenance helper not found: $PROVENANCE_HELPER" >&2; exit 1; }
source "$PROVENANCE_HELPER"

# ---------------------------
# Reference
# ---------------------------
: "${REF_FASTA:?ERROR: REF_FASTA env var not set (export from controller)}"
: "${ERM39_REF_FASTA:?ERROR: ERM39_REF_FASTA env var not set (export from controller)}"

# ---------------------------
# Target gene coordinates (ATCC 19977 reference numbering; 1-based inclusive)
# ---------------------------
RRS_START=1462398; RRS_END=1463901
RRL_START=1464208; RRL_END=1467319
ERM41_START=2345955; ERM41_END=2346476

# erm39 coordinates in M. fortuitum ATCC35855 (CP194182.1)
ERM39_CONTIG="CP194182.1"
ERM39_START=2670342
ERM39_END=2671082
ERM39_POS=2670342
ERM39_LENGTH=741
ERM39_GENE_MIN_DEPTH=10
ERM39_GENE_MIN_PCT=90
ERM39_SITE_MIN_DEPTH=30

# ---------------------------
# Outputs
# ---------------------------
ISO_DIR="${PIPE_OUTDIR}/${ISOLATE}"
MAP_DIR="${ISO_DIR}/mapping"
VAR_DIR="${ISO_DIR}/variants"
CONS_DIR="${VAR_DIR}/consensus"
LOG_DIR="${ISO_DIR}/logs"

ERM39_MAP_DIR="${ISO_DIR}/erm39_mapping"

mkdir -p "$MAP_DIR" "$VAR_DIR" "$CONS_DIR" "$LOG_DIR" "$ERM39_MAP_DIR"
provenance_init "$ISO_DIR"
PROV_DIR="$ISO_DIR/.provenance"
COMMAND_LOG="$PROV_DIR/step2.commands.log"
VERSION_LOG="$PROV_DIR/step2.versions.txt"
: > "$COMMAND_LOG"
: > "$VERSION_LOG"
provenance_section "$COMMAND_LOG" "STEP 2 - Mapping, Variant Calling, Consensus, and erm39 Assessment"

# BWA prints its version in the usage output rather than via --version
BWA_VERSION="$(bwa 2>&1 | awk '/^Version:/{print $2; exit}' || true)"
[[ -n "$BWA_VERSION" ]] || BWA_VERSION="UNAVAILABLE"
printf 'BWA: %s\n' "$BWA_VERSION" >> "$VERSION_LOG"

provenance_version_line "$VERSION_LOG" "SAMtools" samtools --version
provenance_version_line "$VERSION_LOG" "BCFtools" bcftools --version
provenance_version_line "$VERSION_LOG" "Python" "$NTM_PYTHON" --version

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

if [[ ! -s "$ERM39_REF_FASTA" ]]; then
  echo -e "${ISOLATE}\tSTEP2\tFAIL\term39 reference FASTA not found: $ERM39_REF_FASTA" >> "$STATUS_TSV"
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
if [[ ! -s "${ERM39_REF_FASTA}.bwt" ]]; then
  echo -e "${ISOLATE}\tERM39_BWA_INDEX\tFAIL\tMissing ATCC35855 BWA index (.bwt). Controller should index reference first." >> "$STATUS_TSV"
  exit 0
fi
if [[ ! -s "${ERM39_REF_FASTA}.fai" ]]; then
  echo -e "${ISOLATE}\tERM39_FAIDX\tFAIL\tMissing ATCC35855 FASTA index (.fai). Controller should faidx reference first." >> "$STATUS_TSV"
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

provenance_log "$COMMAND_LOG" "bwa mem -t $THREADS $REF_FASTA $R1FASTQ $R2FASTQ 2> ${LOG_DIR}/bwa_mem.${ISOLATE}.err | samtools sort -@ $THREADS -o $BAM - 2> ${LOG_DIR}/samtools_sort.${ISOLATE}.err"

if bwa mem -t "$THREADS" "$REF_FASTA" "$R1FASTQ" "$R2FASTQ" 2> "${LOG_DIR}/bwa_mem.${ISOLATE}.err" \
  | samtools sort -@ "$THREADS" -o "$BAM" - 2> "${LOG_DIR}/samtools_sort.${ISOLATE}.err"; then
  :
else
  echo -e "${ISOLATE}\tMAP\tFAIL\tbwa mem / samtools sort failed (see logs)" >> "$STATUS_TSV"
  exit 0
fi

provenance_log "$COMMAND_LOG" "samtools index $BAM 2> ${LOG_DIR}/samtools_index.${ISOLATE}.err"
samtools index "$BAM" 2> "${LOG_DIR}/samtools_index.${ISOLATE}.err" || {
  echo -e "${ISOLATE}\tMAP\tFAIL\tsamtools index failed (see logs)" >> "$STATUS_TSV"
  exit 0
}

provenance_log "$COMMAND_LOG" "samtools flagstat $BAM > ${MAP_DIR}/flagstat.txt 2> ${LOG_DIR}/flagstat.${ISOLATE}.err"
samtools flagstat "$BAM" > "${MAP_DIR}/flagstat.txt" 2> "${LOG_DIR}/flagstat.${ISOLATE}.err" || true
echo -e "${ISOLATE}\tMAP\tOK\tMapping complete" >> "$STATUS_TSV"

# ---------------------------
# 2b) Variant calling in targets (whole genes)
# ---------------------------
echo "[$(date)] Calling variants in target genes for $ISOLATE" >&2

VCF_GZ="${VAR_DIR}/targets.vcf.gz"

provenance_log "$COMMAND_LOG" "bcftools mpileup -f $REF_FASTA -R $BED -q 20 -Q 20 -a DP,AD -Ou $BAM 2> ${LOG_DIR}/mpileup.${ISOLATE}.err | bcftools call -mv -Oz -o $VCF_GZ 2> ${LOG_DIR}/call.${ISOLATE}.err"

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

provenance_log "$COMMAND_LOG" "bcftools index -t $VCF_GZ 2> ${LOG_DIR}/bcftools_index.${ISOLATE}.err"
bcftools index -t "$VCF_GZ" 2> "${LOG_DIR}/bcftools_index.${ISOLATE}.err" || {
  echo -e "${ISOLATE}\tVCF\tFAIL\tbcftools index failed (see logs)" >> "$STATUS_TSV"
  exit 0
}

DEPTH_TSV="${VAR_DIR}/targets_depth.tsv"
provenance_log "$COMMAND_LOG" "samtools depth -aa -b $BED $BAM > $DEPTH_TSV 2> ${LOG_DIR}/depth.${ISOLATE}.err"
samtools depth -aa -b "$BED" "$BAM" > "$DEPTH_TSV" 2> "${LOG_DIR}/depth.${ISOLATE}.err" || true

RAW_VAR_TSV="${VAR_DIR}/targets_variants.tsv"
provenance_log "$COMMAND_LOG" "bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%AD]\n' $VCF_GZ | awk -v iso='$ISOLATE' -v rrs_s='$RRS_START' -v rrs_e='$RRS_END' -v rrl_s='$RRL_START' -v rrl_e='$RRL_END' -v erm_s='$ERM41_START' -v erm_e='$ERM41_END' 'BEGIN{FS=OFS=\"\t\"; print \"Isolate\",\"Gene\",\"CHROM\",\"POS_ref\",\"POS_gene\",\"REF\",\"ALT\",\"QUAL\",\"DP\",\"AD\"} {chrom=\$1; pos=\$2; ref=\$3; alt=\$4; qual=\$5; dp=\$6; ad=\$7; gene=\"NA\"; gpos=\"NA\"; if (pos>=rrs_s && pos<=rrs_e) {gene=\"rrs\"; gpos=pos-rrs_s+1} else if (pos>=rrl_s && pos<=rrl_e) {gene=\"rrl\"; gpos=pos-rrl_s+1} else if (pos>=erm_s && pos<=erm_e) {gene=\"erm41\"; gpos=pos-erm_s+1} print iso,gene,chrom,pos,gpos,ref,alt,qual,dp,ad}' > $RAW_VAR_TSV"
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
  provenance_log "$COMMAND_LOG" "samtools faidx $REF_FASTA $region 2> ${LOG_DIR}/faidx_region.${ISOLATE}.${gene}.err | bcftools consensus -f - $VCF_GZ 2> ${LOG_DIR}/consensus.${ISOLATE}.${gene}.err | awk -v h='>${ISOLATE}|${gene}|${region}' 'BEGIN{printed=0} /^>/{print h; printed=1; next} {print} END{if(!printed){print h}}' > $out_fa"
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



# ---------------------------
# 2d) erm39 mapping against ATCC35855
#
# BLAST remains the primary erm39 presence detector. This
# mapping analysis evaluates the literature-supported initiation
# codon mutation at erm39 position 1 (G>C; GTG>CTG).
#
# Site evidence is emitted only when >=90% of the 741-bp erm39
# CDS has >=10x coverage against the ATCC35855 reference.
# ---------------------------
echo "[$(date)] Mapping $ISOLATE to ATCC35855 for erm39" >&2

ERM39_BAM="${ERM39_MAP_DIR}/align.sorted.bam"
ERM39_DEPTH="${ERM39_MAP_DIR}/erm39_depth.tsv"
ERM39_PILEUP="${ERM39_MAP_DIR}/erm39_position1.pileup"
ERM39_METRICS="${ERM39_MAP_DIR}/erm39_mapping_metrics.tsv"
ERM39_SITE_TSV="${VAR_DIR}/erm39_site.tsv"

rm -f "$ERM39_SITE_TSV"

provenance_log "$COMMAND_LOG" "bwa mem -t $THREADS $ERM39_REF_FASTA $R1FASTQ $R2FASTQ 2> ${LOG_DIR}/bwa_mem_erm39.${ISOLATE}.err | samtools sort -@ $THREADS -o $ERM39_BAM - 2> ${LOG_DIR}/samtools_sort_erm39.${ISOLATE}.err"

if bwa mem -t "$THREADS" "$ERM39_REF_FASTA" "$R1FASTQ" "$R2FASTQ" \
  2> "${LOG_DIR}/bwa_mem_erm39.${ISOLATE}.err" \
  | samtools sort -@ "$THREADS" -o "$ERM39_BAM" - \
      2> "${LOG_DIR}/samtools_sort_erm39.${ISOLATE}.err"; then
  :
else
  echo -e "${ISOLATE}\tERM39_MAP\tFAIL\tbwa mem / samtools sort failed (see logs)" >> "$STATUS_TSV"
  echo "[$(date)] Step 2 done for $ISOLATE" >&2
  exit 0
fi

provenance_log "$COMMAND_LOG" "samtools index $ERM39_BAM 2> ${LOG_DIR}/samtools_index_erm39.${ISOLATE}.err"
samtools index "$ERM39_BAM" 2> "${LOG_DIR}/samtools_index_erm39.${ISOLATE}.err" || {
  echo -e "${ISOLATE}\tERM39_MAP\tFAIL\tsamtools index failed (see logs)" >> "$STATUS_TSV"
  echo "[$(date)] Step 2 done for $ISOLATE" >&2
  exit 0
}

provenance_log "$COMMAND_LOG" "samtools depth -aa -r ${ERM39_CONTIG}:${ERM39_START}-${ERM39_END} $ERM39_BAM > $ERM39_DEPTH 2> ${LOG_DIR}/depth_erm39.${ISOLATE}.err"
samtools depth \
  -aa \
  -r "${ERM39_CONTIG}:${ERM39_START}-${ERM39_END}" \
  "$ERM39_BAM" \
  > "$ERM39_DEPTH" \
  2> "${LOG_DIR}/depth_erm39.${ISOLATE}.err" \
  || true

provenance_log "$COMMAND_LOG" "samtools mpileup -f $ERM39_REF_FASTA -r ${ERM39_CONTIG}:${ERM39_POS}-${ERM39_POS} $ERM39_BAM > $ERM39_PILEUP 2> ${LOG_DIR}/mpileup_erm39.${ISOLATE}.err"
samtools mpileup \
  -f "$ERM39_REF_FASTA" \
  -r "${ERM39_CONTIG}:${ERM39_POS}-${ERM39_POS}" \
  "$ERM39_BAM" \
  > "$ERM39_PILEUP" \
  2> "${LOG_DIR}/mpileup_erm39.${ISOLATE}.err" \
  || true

provenance_log "$COMMAND_LOG" "$NTM_PYTHON - $ISOLATE $ERM39_DEPTH $ERM39_PILEUP $ERM39_METRICS $ERM39_SITE_TSV $ERM39_LENGTH $ERM39_GENE_MIN_DEPTH $ERM39_GENE_MIN_PCT $ERM39_SITE_MIN_DEPTH <<'PYERM39'"
"$NTM_PYTHON" - \
  "$ISOLATE" \
  "$ERM39_DEPTH" \
  "$ERM39_PILEUP" \
  "$ERM39_METRICS" \
  "$ERM39_SITE_TSV" \
  "$ERM39_LENGTH" \
  "$ERM39_GENE_MIN_DEPTH" \
  "$ERM39_GENE_MIN_PCT" \
  "$ERM39_SITE_MIN_DEPTH" <<'PYERM39'
import sys

(
    isolate,
    depth_path,
    pileup_path,
    metrics_path,
    site_path,
    erm39_length,
    gene_min_depth,
    gene_min_pct,
    site_min_depth,
) = sys.argv[1:]

erm39_length = int(erm39_length)
gene_min_depth = int(gene_min_depth)
gene_min_pct = float(gene_min_pct)
site_min_depth = int(site_min_depth)

ref_base = "G"
alt_base = "C"

depths = []
try:
    with open(depth_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) >= 3:
                depths.append(int(fields[2]))
except FileNotFoundError:
    pass

if len(depths) < erm39_length:
    depths.extend([0] * (erm39_length - len(depths)))
elif len(depths) > erm39_length:
    depths = depths[:erm39_length]

covered_gene = sum(d >= gene_min_depth for d in depths)
pct_covered = 100.0 * covered_gene / erm39_length
gene_callable = pct_covered >= gene_min_pct

counts = {"A": 0, "C": 0, "G": 0, "T": 0}

try:
    with open(pileup_path) as handle:
        line = handle.readline().rstrip("\n")
except FileNotFoundError:
    line = ""

if line:
    fields = line.split("\t")
    if len(fields) >= 5:
        read_bases = fields[4]
        i = 0
        while i < len(read_bases):
            c = read_bases[i]
            if c == "^":
                i += 2
                continue
            if c == "$":
                i += 1
                continue
            if c in ".,":
                counts[ref_base] += 1
                i += 1
                continue
            if c.upper() in counts:
                counts[c.upper()] += 1
                i += 1
                continue
            if c in "+-":
                i += 1
                number = ""
                while i < len(read_bases) and read_bases[i].isdigit():
                    number += read_bases[i]
                    i += 1
                if number:
                    i += int(number)
                continue
            i += 1

depth = sum(counts.values())
c_af = counts[alt_base] / depth if depth > 0 else 0.0

if depth < site_min_depth:
    site_call = "INDETERMINATE"
elif c_af >= 0.90:
    site_call = "MUT"
elif c_af >= 0.10:
    site_call = "MIXED"
else:
    site_call = "WT"

with open(metrics_path, "w") as out:
    out.write(
        "Isolate\tGene\tgene_length\tgene_min_depth\t"
        "pct_gene_covered\tgene_callable\tposition\t"
        "ref_position\tREF\tALT\tDepth\tA\tC\tG\tT\tAF\tsite_call\n"
    )
    out.write(
        f"{isolate}\term39\t{erm39_length}\t{gene_min_depth}\t"
        f"{pct_covered:.2f}\t{str(gene_callable).upper()}\t1\t"
        f"2670342\t{ref_base}\t{alt_base}\t{depth}\t"
        f"{counts['A']}\t{counts['C']}\t{counts['G']}\t{counts['T']}\t"
        f"{c_af:.4f}\t{site_call}\n"
    )

if gene_callable:
    with open(site_path, "w") as out:
        out.write(
            "Isolate\tGene\tposition\tDepth\tREF\tALT\tQUAL\tDP\tAD\tAF\n"
        )
        out.write(
            f"{isolate}\term39\t1\t{depth}\t{ref_base}\t{alt_base}\t\t"
            f"{depth}\t{counts[ref_base]},{counts[alt_base]}\t{c_af:.4f}\n"
        )
PYERM39

echo -e "${ISOLATE}\tERM39_MAP\tOK\term39 mapping assessment complete" >> "$STATUS_TSV"

echo "[$(date)] Step 2 done for $ISOLATE" >&2


