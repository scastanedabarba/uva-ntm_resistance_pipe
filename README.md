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

    git clone https://github.com/<scastanedabarba>/uva-ntm_resistance_pipe.git
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

```{=html}
<!-- -->
```
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

-   rrs: 1,462,398 -- 1,463,901\
-   rrl: 1,464,208 -- 1,467,319\
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

summary/blast_top_hits.tsv\
summary/sites_evidence.tsv\
summary/erm41_truncation_metrics.tsv\
summary/interpretation.tsv\
summary/myco_prediction_summary.xlsx

The Excel file aggregates variant, truncation, and BLAST summaries for
review but TSV files remain the primary reproducible outputs.

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
