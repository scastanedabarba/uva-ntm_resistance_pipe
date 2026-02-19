# UVA Mycobacterial Resistance Pipeline

## Overview

This pipeline performs:

1. Assembly + BLAST detection of erm genes  
2. Targeted variant calling (rrl, rrs, erm41)  
3. Coverage-based erm41 truncation assessment  
4. Interpretation of clarithromycin and amikacin resistance  

It supports:
- Clinical isolates
- Simulated ATCC19977 datasets

---

# ⚠️ IMPORTANT: Simulation vs Pipeline

## Simulation DOES NOT run the pipeline

Running the simulation scripts **only generates synthetic FASTQ datasets**.

It does NOT:
- assemble reads
- call variants
- generate resistance predictions

After simulation, you must run the full pipeline on the simulated dataset.

---

# Part 1 – Generate Simulated Dataset (Optional)

Example:

```bash
bash simulation/run_atcc_sim_pipeline.sh \
  --outdir /scratch/sgj4qr/atcc_sim_test \
  --coverage 100
```

This generates:

```
/scratch/sgj4qr/atcc_sim_test/
    isolate1/
        isolate1_R1.fastq.gz
        isolate1_R2.fastq.gz
    isolate2/
    ...
```

This step only creates FASTQs.

---

# Part 2 – Run the Pipeline

## Clinical dataset

```bash
bash submit_ntm_pipeline.sh \
  linelist.tsv \
  /scratch/sgj4qr/myco_run
```

## Simulated dataset

```bash
bash submit_ntm_pipeline.sh \
  linelist.tsv \
  /scratch/sgj4qr/atcc_sim_test
```

Where:
- `linelist.tsv` contains isolate IDs (first column only)
- `outdir` is the root pipeline directory

---

# Pipeline Structure

For each isolate:

```
outdir/
  ISOLATE/
    assembly/
    blast/
    mapping/
    variants/
    status_step1.tsv
    status_step2.tsv
```

---

# Output Files Explained

All final outputs are written to:

```
outdir/summary/
```

---

## 1️⃣ blast_top_hits.tsv

Summarizes top BLAST hit per gene group:

- erm41
- erm39
- erm55

Columns:
- Pident
- QcovPct
- TotalHits
- PreferredHits (≥90% identity AND ≥90% coverage)

Used for:
- Detecting erm39 / erm55
- erm41 presence (fallback if coverage insufficient)

---

## 2️⃣ sites_evidence.tsv

Contains only variant sites observed at positions of interest:

| Gene  | Position     | Meaning |
|--------|-------------|---------|
| rrl    | 2269–2293   | Macrolide resistance mutations |
| rrs    | 1373–1458   | Amikacin resistance mutations |
| erm41  | 19, 28      | Inducible macrolide resistance genotype |

Columns:
- Depth (reference coordinate depth)
- REF
- ALT
- DP
- AD
- AF (allele frequency)

Allele Frequency logic:
- ≥0.90 → MUT
- 0.10–0.89 → MIXED
- <0.10 → WT

---

## 3️⃣ erm41_truncation_metrics.tsv

Coverage-based truncation assessment.

Key column:
- `del_to_flank_ratio`

Logic:
If flank median depth ≥ threshold AND  
   deletion-region median depth / flank median depth ≤ 0.10  
→ TRUNCATED  
Else → NOT_TRUNCATED  

If flanks not callable → INDETERMINATE

---

## 4️⃣ interpretation.tsv  ⭐ Final Calls

Columns include:

- clarithromycin_call
- clarithromycin_reason
- amikacin_call
- amikacin_reason
- Individual site states
- erm41_truncation
- erm39_detected
- erm55_detected
- notes

This is the authoritative final output.

---

# Interpretation Logic

## Clarithromycin

Priority order:

1. Any rrl MUT → Resistant  
2. rrl MIXED → Resistance possible  
3. rrl WT → evaluate erm41  

If erm41:

- TRUNCATED → Susceptible  
- NOT_TRUNCATED:
  - 28 = C → Susceptible  
  - 28 = T & 19 = C → Resistant  
  - Mixed → Resistance possible  

erm39 / erm55 detected:
→ Forces "Resistance possible" unless already Resistant

---

## Amikacin

- Any rrs MUT → Resistant  
- Mixed → Resistance possible  
- Low depth → Indeterminate  
- Otherwise → Susceptible  

---

# Excel Summary

```
summary/myco_prediction_summary.xlsx
```

Tabs:
- variants
- truncation
- blast

Convenience file only — interpretation.tsv is authoritative.

---

# Status Tracking

Each isolate contains:

## status_step1.tsv
Assembly + BLAST

## status_step2.tsv
Mapping + variant calling

---

# Rerunning Failed Isolates

If an isolate failed:

```bash
bash submit_ntm_pipeline.sh \
  linelist.tsv \
  /scratch/sgj4qr/myco_run \
  --isolate ISOLATE_ID
```

This reruns only that isolate.

---

# Rerun Only Summary + Interpretation

After fixing failed isolates:

```bash
sbatch \
  --export=ALL,REPO_ROOT=/scratch/sgj4qr/mycobac_validation/scripts \
  workflow/03_compile_and_interpret.slurm \
  /scratch/sgj4qr/myco_run \
  linelist.tsv \
  references/nucleotide.fna \
  CU458896.1
```

This regenerates:

- blast_top_hits.tsv
- sites_evidence.tsv
- erm41_truncation_metrics.tsv
- interpretation.tsv
- myco_prediction_summary.xlsx
- combined status tables

No re-assembly occurs.

