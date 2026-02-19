
# UVA NTM Resistance Pipeline

## Purpose

**uva-ntm_resistance_pipe** is a reproducible SLURM-based pipeline for detecting macrolide
and aminoglycoside resistance in *Mycobacterium abscessus* using whole genome sequencing data.

The pipeline performs:

1. Assembly (SPAdes)  
2. Target gene BLAST detection  
3. Reference mapping (ATCC 19977)  
4. Variant calling at resistance-associated loci  
5. erm41 truncation detection (coverage-based)  
6. Structured drug resistance interpretation  

Primary outputs include variant evidence tables, gene detection summaries, truncation metrics, and final antibiotic interpretations.

---

# Installation

## 1) Clone Repository

```bash
git clone https://github.com/scastanedabarba/uva-ntm_resistance_pipe.git
cd mycobac_validation
```

## 2) Required Software

Tested on UVA Rivanna HPC.

- SPAdes ≥ 4.2.0  
- bwa ≥ 0.7.17  
- samtools ≥ 1.15  
- bcftools ≥ 1.15  
- BLAST+  
- Python ≥ 3.9  

### Python Packages

```bash
pip install openpyxl matplotlib
```

---

# Clinical Dataset (Primary Use Case)

## Required Inputs

- Linelist (CSV or TSV)
- Trimmed paired-end reads located at:

```
/project/amr_services/qc/<run>/<isolate>/<isolate>_R1.trim.fq.gz
/project/amr_services/qc/<run>/<isolate>/<isolate>_R2.trim.fq.gz
```

## Example Linelist (Two Column Format)

```text
isolate	run
VALID_0001	251216_M70741_0293_000000000-M84N5
VALID_0002	251216_M70741_0293_000000000-M84N5
```

## Run Pipeline

```bash
bash bin/submit_ntm_pipeline.sh linelist.tsv /path/to/output
```

---

# Reference Genome and Coordinate System

All mapping is performed against:

**Mycobacterium abscessus ATCC 19977**  
GenBank: CU458896.1  
Associated assembly: OQ656457.1  

## Gene Coordinates (Genome-Based, 1-Based)

- rrs: 1,462,398–1,463,901  
- rrl: 1,464,208–1,467,319  
- erm41: 2,345,955–2,346,476  

## Site-of-Interest Positions (Gene-Relative, 1-Based)

These positions are **1-based relative to the gene start**, not whole-genome positions.

- **rrl:** 2269, 2270, 2271, 2281, 2293  
- **rrs:** 1373, 1375, 1376, 1458  
- **erm41:** 19, 28  

Conversion formula:

```
reference_position = gene_start + (gene_position - 1)
```

This removes ambiguity between genome-based and gene-relative numbering.

---

# Interpretation Logic

## Clarithromycin

1. Any rrl MUT → Resistant  
2. rrl MIXED → Resistance possible  
3. rrl WT → evaluate erm41:
   - Truncated → Susceptible  
   - Full-length + 28=C → Susceptible  
   - Full-length + 28=T + 19=C → Resistant  
   - Mixed at 28 or 19 → Resistance possible  

## Amikacin

1. rrs MUT → Resistant  
2. rrs MIXED → Resistance possible  
3. Insufficient depth → Indeterminate  
4. All WT → Susceptible  

---

# Output Files

All final outputs are written to:

```
<outdir>/summary/
```

The TSV files are the primary reproducible outputs.  
The Excel file is a convenience review file.

---

## 1. summary/blast_top_hits.tsv

**Purpose:**  
Summarizes detection of erm41, erm39, and erm55 genes using BLAST.

**Detection Threshold:**  
Percent identity ≥ 90% AND Query coverage ≥ 90%

| Column | Description |
|--------|-------------|
| Isolate | Isolate ID |
| GeneGroup | erm41 / erm39 / erm55 |
| QueryID | Query FASTA sequence ID |
| QueryLen | Length of query gene |
| Subject | Contig hit in assembly |
| Pident | Percent identity |
| AlnLen | Alignment length |
| QcovPct | Percent query coverage |
| Qstart, Qend | Query alignment coordinates |
| Sstart, Send | Subject alignment coordinates |
| Evalue | BLAST E-value |
| Bitscore | BLAST bitscore |
| TotalHits | Total hits detected |
| PreferredHits | Hits ≥90% identity AND ≥90% coverage |

---

## 2. summary/sites_evidence.tsv

**Purpose:**  
Variant evidence at predefined resistance-associated positions.

| Column | Description |
|--------|-------------|
| Isolate | Isolate ID |
| Gene | rrl / rrs / erm41 |
| position | Gene-relative 1-based position |
| Depth | Read depth |
| REF | Reference base |
| ALT | Alternate base |
| QUAL | Variant quality |
| DP | Total depth |
| AD | Allele depths |
| AF | Allele frequency |

**Allele Frequency Interpretation:**

- AF ≥ 0.90 → MUT  
- 0.10 ≤ AF < 0.90 → MIXED  
- AF < 0.10 → WT  
- Depth < 30× → INDETERMINATE  

---

## 3. summary/erm41_truncation_metrics.tsv

**Purpose:**  
Detect erm41 truncation via coverage analysis.

Deletion inferred if:

```
Median depth (deletion region) / Median depth (flanks) ≤ 0.10
AND flank coverage ≥ 20×
```

| Column | Description |
|--------|-------------|
| Isolate | Isolate ID |
| erm_len | Gene length |
| erm_median_depth | Median gene depth |
| left_median_depth | Median left flank |
| right_median_depth | Median right flank |
| flank_median_depth | Median flank depth |
| erm_del_range | Gene-relative deletion region |
| del_median_depth | Deletion region depth |
| del_to_flank_ratio | Deletion depth ÷ flank depth |
| erm41_callable | TRUE if flank depth ≥ 20× |

---

## 4. summary/interpretation.tsv

**Purpose:**  
Final phenotype prediction for clarithromycin and amikacin.

| Column | Description |
|--------|-------------|
| Isolate | Isolate ID |
| clarithromycin_call | Final clarithromycin interpretation |
| clarithromycin_reason | Rule applied |
| amikacin_call | Final amikacin interpretation |
| amikacin_reason | Rule applied |
| rrl_* | Site-level calls |
| erm41_28 | Base at position 28 |
| erm41_19 | Base at position 19 |
| rrs_* | Site-level calls |
| erm41_truncation | TRUE/FALSE |
| erm39_detected | TRUE/FALSE |
| erm55_detected | TRUE/FALSE |
| notes | Additional logic flags |

---

## 5. summary/myco_prediction_summary.xlsx

Convenience aggregation of variant, truncation, and BLAST results.

The TSV files remain the authoritative reproducible outputs.

---

# Simulated ATCC Dataset

## Purpose

Validation and benchmarking dataset containing:

- WT reference  
- Point mutants  
- Deletion mutants  
- Mixed allele datasets  
- High and low coverage reads  

## Step 1 – Generate Dataset

```bash
bash bin/run_atcc_simulation.sh --outdir /path/to/output
```

Creates:

```
/path/to/output/atcc_dataset/
```

This step generates data only.

## Step 2 – Run Pipeline on Simulation

```bash
bash bin/submit_ntm_pipeline.sh   /path/to/output/atcc_dataset/linelist.tsv   /path/to/output   --simulated
```

---

# References

ATCC 19977 chromosome: CU458896.1, OQ656457.1  
erm55: OQ656455.1, OQ656456.1, OQ656457.1  
erm39: AY487229.1  

Key literature:  
https://www.sciencedirect.com/science/article/pii/S0167701217302749  
https://www.sciencedirect.com/science/article/pii/S1525157821002592  

Software used: SPAdes, BWA, SAMtools, BCFtools, BLAST+, Python

