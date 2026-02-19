# UVA NTM Resistance Pipeline

## Purpose

uva-ntm_resistance_pipe is a reproducible SLURM-based pipeline for detecting macrolide
and aminoglycoside resistance in *Mycobacterium abscessus* using whole
genome sequencing data.

The pipeline performs:

1.  Assembly (SPAdes)
2.  Target gene BLAST detection
3.  Reference mapping (ATCC 19977)
4.  Variant calling at resistance-associated loci
5.  erm41 truncation detection (coverage-based)
6.  Structured drug resistance interpretation

Primary outputs include variant evidence tables, gene detection
summaries, truncation metrics, and final antibiotic interpretations.

------------------------------------------------------------------------

# Installation

## 1) Clone Repository

    git clone https://github.com/scastanedabarba/uva-ntm_resistance_pipe.git
    cd mycobac_validation

## 2) Required Software

Tested on UVA Rivanna.

-   SPAdes ≥ 4.2.0\
-   bwa ≥ 0.7.17\
-   samtools ≥ 1.15\
-   bcftools ≥ 1.15\
-   BLAST+\
-   Python ≥ 3.9

### Python Packages

-   openpyxl\
-   matplotlib

------------------------------------------------------------------------

# Clinical Dataset (Primary Use Case)

## Required Inputs

-   Linelist (CSV or TSV)
-   Trimmed paired-end reads located at:

    /project/amr_services/qc/<run>/<isolate>/<isolate>_R1.trim.fq.gz
    /project/amr_services/qc/<run>/<isolate>/<isolate>_R2.trim.fq.gz

## Example Linelist

Two-column format:

    isolate run
    VALID_0001  251216_M70741_0293_000000000-M84N5
    VALID_0002  251216_M70741_0293_000000000-M84N5

## Run Pipeline

    bash bin/submit_ntm_pipeline.sh linelist.tsv /path/to/output

------------------------------------------------------------------------

# Reference Genome and Coordinate System

All mapping is performed against:

**Mycobacterium abscessus ATCC 19977**\
GenBank: CU458896.1\
Associated assembly: OQ656457.1

## Gene Coordinates (1-based reference positions)

-   rrs: 1,462,398 -- 1,463,901
-   rrl: 1,464,208 -- 1,467,319
-   erm41: 2,345,955 -- 2,346,476

## Site-of-Interest Positions (Gene-Relative)

These are **1-based positions relative to the start of each gene**, not
genome coordinates:

-   rrl: 2269, 2270, 2271, 2281, 2293
-   rrs: 1373, 1375, 1376, 1458
-   erm41: 19, 28

Internal conversion:

    reference_position = gene_start + (gene_position - 1)

This ensures reproducibility and removes ambiguity between gene-relative
and genome-relative numbering.

------------------------------------------------------------------------

# Interpretation Logic

## Clarithromycin

1.  rrl mutation → Resistant
2.  rrl mixed → Resistance possible
3.  rrl WT → evaluate erm41:
    -   Truncated → Susceptible
    -   Full-length + 28=C → Susceptible
    -   Full-length + 28=T + 19=C → Resistant
    -   Mixed at 28 or 19 → Resistance possible

## Amikacin

1.  rrs mutation → Resistant
2.  rrs mixed → Resistance possible
3.  Insufficient depth → Indeterminate
4.  All WT → Susceptible

------------------------------------------------------------------------

# Output Files

All final outputs are written to:

<outdir>/summary/

The TSV files are the primary reproducible outputs.
The Excel file is a convenience report.

---------------------------------------------------------------------

## 1. summary/blast_top_hits.tsv

Purpose:
Summarizes detection of erm41, erm39, and erm55 genes using BLAST
against assembled contigs.

Used to:
- Confirm gene presence
- Support erm41 applicability logic
- Flag potential macrolide resistance from erm39/erm55

Columns:
Isolate              : Isolate ID
GeneGroup            : erm41 / erm39 / erm55
QueryID              : Query FASTA sequence ID
QueryLen             : Length of query gene
Subject              : Contig hit in assembly
Pident               : Percent identity
AlnLen               : Alignment length
QcovPct              : Percent query coverage
Qstart, Qend         : Query alignment coordinates
Sstart, Send         : Subject alignment coordinates
Evalue               : BLAST E-value
Bitscore             : BLAST bitscore
TotalHits            : Total hits detected
PreferredHits        : Hits meeting ≥90% identity AND ≥90% coverage

Detection Threshold:
A gene is considered detected if:
Percent identity ≥ 90%
AND
Query coverage ≥ 90%

---------------------------------------------------------------------

## 2. summary/sites_evidence.tsv

Purpose:
Lists variant evidence at predefined resistance-associated positions.

Important:
Positions are 1-based relative to the gene start,
NOT relative to the whole ATCC genome.

Sites Evaluated:

rrl:   2269, 2270, 2271, 2281, 2293
rrs:   1373, 1375, 1376, 1458
erm41: 19, 28

Columns:
Isolate     : Isolate ID
Gene        : rrl / rrs / erm41
position    : Gene-relative 1-based position
Depth       : Read depth at corresponding reference coordinate
REF         : Reference base
ALT         : Alternate base
QUAL        : Variant quality
DP          : Total depth from VCF
AD          : Allele depths (ref,alt)
AF          : Allele frequency (alt_count / depth)

Allele Frequency Interpretation:

AF ≥ 0.90         : MUT
0.10 ≤ AF < 0.90  : MIXED
AF < 0.10         : WT
Depth < threshold : INDETERMINATE

Default depth threshold: 30×.

---------------------------------------------------------------------

## 3. summary/erm41_truncation_metrics.tsv

Purpose:
Determines whether erm41 is truncated based on coverage.

Strategy:
Coverage is evaluated across:

Left flank:      gene positions 20–140
Right flank:     gene positions 450–560
Deletion region: gene positions 159–432

Deletion inferred if:

Median depth (deletion region)
/
Median depth (flanks)
≤ 0.10

AND flanking coverage ≥ 20×.

Columns:
Isolate               : Isolate ID
erm_len               : Length of erm41 gene
erm_median_depth      : Median depth across full gene
left_median_depth     : Median depth left flank
right_median_depth    : Median depth right flank
flank_median_depth    : Median of left + right
erm_del_range         : Gene-relative deletion region
del_median_depth      : Median depth in deletion region
del_to_flank_ratio    : Deletion depth ÷ flank depth
erm41_callable        : TRUE if flank depth ≥ 20×

---------------------------------------------------------------------

## 4. summary/interpretation.tsv

Purpose:
Final phenotype prediction per isolate for:

- Clarithromycin
- Amikacin

Integrates:
- rrl mutations
- rrs mutations
- erm41 truncation
- erm41 genotype (positions 19, 28)
- erm39 detection
- erm55 detection

Columns:

Isolate
clarithromycin_call
clarithromycin_reason
amikacin_call
amikacin_reason
rrl_* site calls
erm41_28
erm41_19
rrs_* site calls
erm41_truncation
erm39_detected
erm55_detected
notes

Clarithromycin Logic (Simplified):

1. Any rrl MUT → Resistant
2. rrl MIXED → Resistance possible
3. rrl WT → evaluate erm41
   - Truncated → Susceptible
   - Full-length + 28=T & 19=C → Resistant
   - Full-length + 28=T & 19=T → Susceptible
   - Mixed → Resistance possible
4. erm39/erm55 detected → may upgrade to Resistance possible

Amikacin Logic:

1. rrs MUT → Resistant
2. rrs MIXED → Resistance possible
3. Insufficient depth → Indeterminate
4. Otherwise → Susceptible

---------------------------------------------------------------------

## 5. summary/myco_prediction_summary.xlsx

Purpose:
Convenience aggregation of:

- Variant evidence
- Truncation metrics
- BLAST summary

Intended for clinician-friendly review.

The TSV files remain the primary reproducible outputs.

---------------------------------------------------------------------

Reproducibility Note:

All phenotype decisions are derived exclusively from:

- sites_evidence.tsv
- erm41_truncation_metrics.tsv
- blast_top_hits.tsv

interpretation.tsv is a deterministic transformation of these files.

------------------------------------------------------------------------

# Simulated ATCC Dataset

## Purpose

The simulation dataset is used for validation and benchmarking. It
generates:

-   WT reference
-   Point mutants
-   Deletion mutants
-   Mixed allele datasets
-   High and low coverage reads

## Step 1 -- Generate Dataset

    bash bin/run_atcc_simulation.sh --outdir /path/to/output

This creates:

    /path/to/output/atcc_dataset/

This step **only generates data**. It does not run resistance analysis.

## Step 2 -- Run Pipeline on Simulation

    bash bin/submit_ntm_pipeline.sh \
      /path/to/output/atcc_dataset/linelist.tsv \
      /path/to/output \
      --simulated

The --simulated flag instructs the pipeline to use:

    <outdir>/atcc_dataset/reads/<isolate>/

instead of the clinical QC path.

------------------------------------------------------------------------

# References

ATCC 19977 chromosome: CU458896.1, OQ656457.1\
erm55: OQ656455.1, OQ656456.1, OQ656457.1\
erm39: AY487229.1

Key literature: -
https://www.sciencedirect.com/science/article/pii/S0167701217302749 -
https://www.sciencedirect.com/science/article/pii/S1525157821002592

Software used: - SPAdes - BWA - SAMtools - BCFtools - BLAST+ - Python
(Biopython, openpyxl)
