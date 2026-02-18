# UVA NTM Resistance Pipeline

Pipeline for genomic prediction of **clarithromycin and amikacin resistance** in *Mycobacterium abscessus complex* using whole-genome sequencing data.

Designed for SLURM-based HPC environments (e.g., UVA Rivanna).

---

## Overview

This pipeline performs:

1. Assembly-based detection of macrolide resistance genes (erm41, erm39, erm55)
2. Reference-based variant calling of:
   - rrl
   - rrs
   - erm41
3. Coverage-based detection of erm41 truncation
4. Rule-based antimicrobial susceptibility interpretation
5. Optional ATCC 19977 simulation framework for validation

---

## Repository Structure

bin/ Entry-point wrapper scripts
workflow/ SLURM step scripts (01–03)
scripts/ Python logic (compile + interpret)
simulation/ ATCC validation dataset generator
references/ Packaged reference FASTAs


---

## Requirements

### System
- Linux
- SLURM scheduler

### Required software (via modules or conda)
- spades
- bwa
- samtools
- bcftools
- blastn
- wgsim
- python3 (with biopython, openpyxl)

---

## References Included

references/ATCC19977.fasta
references/nucleotide.fna


- ATCC 19977 reference genome
- Target gene FASTA for BLAST (erm variants)

Indexes (bwa, samtools) are generated at runtime and are NOT stored in the repo.

---

## Running the Pipeline (Real Data)

bash bin/submit_ntm_pipeline.sh
<linelist.tsv>
<workdir>
<outdir>


### Linelist formats

- 2-column: isolate, run
- OR 1-column isolate (requires --reads-root)

### Output structure

Per isolate:

assembly/
blast/
mapping/
variants/
logs/


Global outputs:

summary/
blast_top_hits.tsv
sites_evidence.tsv
erm41_truncation_metrics.tsv
interpretation.tsv
myco_prediction_summary.xlsx


---

## Running Simulation Validation

Generate mutant FASTAs, simulate reads, and test pipeline logic:

bash bin/run_atcc_simulation.sh --workdir <dir>


This will:

1. Generate SNP and deletion mutants
2. Simulate high-coverage reads
3. Simulate low-coverage datasets
4. Create mixed allele datasets (5% and 50%)

---

## Interpretation Logic

Clarithromycin decision tree prioritizes:

1. rrl mutations
2. erm41 truncation
3. erm41 19/28 genotype
4. erm39/erm55 detection

Amikacin decision based on rrs mutations.

Default depth threshold: 30x (configurable).

---

## Versioning

Tag stable releases:

git tag -a v0.1.0 -m "Initial stable release"
git push origin v0.1.0


---

## Author

Salvador Castaneda  
University of Virginia  
AMR Services

