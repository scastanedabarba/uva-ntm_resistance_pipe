# UVA NTM Resistance Pipeline

## Purpose

**uva-ntm_resistance_pipe** is a reproducible SLURM-based whole-genome sequencing (WGS)
pipeline for detecting and interpreting genetic determinants associated with
clarithromycin and amikacin resistance in rapidly growing nontuberculous
mycobacteria (NTM).

The pipeline has been developed and validated for resistance-associated loci
relevant to *Mycobacterium abscessus*, *Mycobacterium fortuitum*, and
*Mycobacterium chelonae*.

The pipeline performs:

1. De novo assembly using SPAdes
2. BLAST detection of `erm41`, `erm39`, and `erm55`
3. Reference mapping for resistance-associated loci
4. Variant calling at predefined `rrl`, `rrs`, and `erm41` positions
5. Coverage-based detection of `erm41` truncation
6. Targeted assessment of the `erm39` initiation codon
7. Structured clarithromycin and amikacin resistance interpretation
8. Generation of reproducible TSV outputs and an Excel evidence summary with target-site coverage

The primary outputs are:

- `interpretation.tsv` — final WGS-based clarithromycin and amikacin interpretations
- `myco_prediction_summary.xlsx` — supporting evidence used to review the final interpretations

Additional TSV files containing BLAST, site-level coverage/variant evidence, and
truncation evidence are retained for reproducibility and detailed review.

---

# Installation

## 1. Clone Repository

```bash
git clone https://github.com/scastanedabarba/uva-ntm_resistance_pipe.git
```

## 2. Software Environment

The pipeline is designed for execution on UVA Rivanna HPC.

Required command-line software includes:

- SPAdes 4.2.0
- BWA 0.7.17
- SAMtools 1.23
- BCFtools 1.23
- BLAST+ 2.11
- Python 3.11

Python dependencies include:

- Biopython
- NumPy
- pandas
- matplotlib
- openpyxl

Software environments are managed through:

```text
config/setup_environment.sh
```

The pipeline uses the dedicated Python environment:

```text
/project/amr_services/.conda/ntm_resistance_pipe
```

and project software modules for compiled command-line tools.

---

# Pipeline Execution

## Required Inputs

The primary input is a CSV or TSV linelist.

The first column must contain the isolate identifier. For routine AMR Services
data, the sequencing run must also be provided.

Example:

```text
isolate,run
VALID_0001,251216_M70741_0293_000000000-M84N5
VALID_0002,251216_M70741_0293_000000000-M84N5
```

Trimmed paired-end reads are expected at:

```text
/project/amr_services/qc/<run>/<isolate>/<isolate>_R1.trim.fq.gz
/project/amr_services/qc/<run>/<isolate>/<isolate>_R2.trim.fq.gz
```

## Run Pipeline

```bash
bash bin/submit_ntm_pipeline.sh \
    linelist.tsv \
    /path/to/output
```

Help and available options can be displayed with:

```bash
bash bin/submit_ntm_pipeline.sh --help
```

Available options include:

```text
--simulated <DIR>  Run in ATCC simulated mode using reads in:
                   <DIR>/atcc_dataset/reads/<isolate>/<isolate>_R{1,2}.trim.fq.gz

--isolate <ID>     Run only a single isolate
                   (must be present in linelist)

--partition <p>    SLURM partition
                   (default: standard)

--ref-fasta <p>    Override the primary reference FASTA
                   (default: repo/references/ATCC19977.fasta)

-h, --help         Show help and exit
```

The submission script performs pipeline setup and submits the individual
SLURM workflow jobs with the required dependencies.

---

# Pipeline Workflow

For each isolate, the pipeline performs three major steps.

## Step 1: Assembly and Resistance Gene Detection

Trimmed paired-end reads are assembled using SPAdes.

The resulting assembly is searched using BLAST for:

- `erm41`
- `erm39`
- `erm55`

BLAST detection is used to determine the presence of these
resistance-associated genes.

## Step 2: Reference Mapping and Resistance Locus Assessment

Reads are mapped to reference genomes to evaluate predefined
resistance-associated loci.

### ATCC 19977 Mapping

Reads are mapped against the *M. abscessus* ATCC 19977 reference for
assessment of:

- `rrl`
- `rrs`
- `erm41`

This mapping is used for:

- site-specific variant detection
- sequencing-depth assessment
- `erm41` genotype determination
- `erm41` truncation assessment

### ATCC 35855 Mapping

Reads are additionally mapped against the *M. fortuitum* ATCC 35855 reference
for targeted assessment of `erm39`.

BLAST remains the primary method for determining `erm39` presence. Reference
mapping is used to assess the characterized susceptibility-associated
initiation-codon variant when sufficient `erm39` mapping coverage is available.

## Step 3: Compilation and Interpretation

Evidence from assembly, BLAST, mapping, variant detection, and coverage
assessment is compiled into final clarithromycin and amikacin interpretations.

The final interpretation table and Excel evidence workbook are generated in
this step.

---

# Reference Genomes and Coordinate Systems

## Mycobacterium abscessus ATCC 19977

Primary reference:

```text
CU458896.1
```

### Gene Coordinates

Genome coordinates are 1-based.

| Gene | Reference coordinates |
|---|---:|
| `rrs` | 1,462,398–1,463,901 |
| `rrl` | 1,464,208–1,467,319 |
| `erm41` | 2,345,955–2,346,476 |

### Resistance-Associated Positions

Positions reported by the pipeline are 1-based relative to the corresponding
gene.

| Gene | Positions evaluated |
|---|---|
| `rrl` | 2269, 2270, 2271, 2281, 2293 |
| `rrs` | 1373, 1375, 1376, 1458 |
| `erm41` | 19, 28 |

For a forward-strand gene:

```text
reference_position = gene_start + (gene_position - 1)
```

This distinction is important because reported resistance positions use
gene-relative numbering rather than whole-genome coordinates.

## Mycobacterium fortuitum ATCC 35855

The complete ATCC 35855 reference is used for targeted `erm39` mapping.

`erm39` is located at:

```text
Reference contig: CP194182.1
Coordinates:      2,670,342–2,671,082
Length:           741 bp
Strand:           forward
```

The `erm39` initiation codon is:

```text
Gene positions:       1–3
Reference sequence:   GTG
Reference position 1: CP194182.1:2670342
```

The pipeline specifically evaluates the G→C substitution at gene position 1,
corresponding to the characterized GTG→CTG initiation-codon change.

---

# Interpretation Logic

The pipeline assigns WGS-based interpretations using predefined resistance
determinants and variants included within the validated rule set.

## Allele Frequency Thresholds

For callable variant positions:

```text
AF < 0.10          Reference
0.10 ≤ AF < 0.90  Mixed
AF ≥ 0.90          Mutant
```

Target positions require:

```text
Depth ≥ 30×
```

Sites below the required depth are reported as indeterminate rather than being
assumed to match the reference.

---

## Clarithromycin

Clarithromycin interpretation incorporates evidence from:

- `rrl`
- `erm41`
- `erm39`
- `erm55`

### rrl

Any mutant allele at a predefined `rrl` resistance-associated position is
interpreted as:

```text
Resistant
```

A mixed allele can result in:

```text
Resistance possible
```

depending on the remaining resistance evidence.

### erm41

When `erm41` is applicable and sufficiently covered:

- Truncated `erm41` → Susceptible in the absence of an overriding resistance determinant
- Full-length `erm41`, position 28 = C → Susceptible
- Full-length `erm41`, position 28 = T and position 19 = C → Resistant
- Full-length `erm41`, position 28 = T and position 19 = T → Susceptible
- Mixed relevant alleles → Resistance possible
- Insufficient evidence → Indeterminate when an interpretation cannot otherwise be made

### erm39

`erm39` is first detected using BLAST.

When `erm39` is detected, the pipeline evaluates its initiation codon when
mapping is sufficiently callable.

The reference initiation codon is:

```text
GTG
```

The characterized susceptibility-associated mutation is:

```text
GTG → CTG
```

Interpretation includes:

- Intact GTG initiation codon → Resistant
- CTG initiation-codon mutant → Susceptible in the absence of an overriding resistance determinant
- Mixed G/C population → Resistance possible
- `erm39` detected but initiation codon not callable → Resistant by default

The pipeline does not attempt to assign phenotypes to all possible sequence
variation within `erm39`. Only variants incorporated into the predefined
interpretation rules are used for phenotype assignment.

### erm55

Detection of `erm55` provides additional evidence of macrolide resistance and
can result in a `Resistance possible` interpretation when a definitive
resistant interpretation has not already been assigned.

---

## Amikacin

Amikacin interpretation is based on predefined resistance-associated `rrs`
positions.

- Any mutant allele at a target `rrs` position → Resistant
- Mixed allele at a target position → Resistance possible
- Insufficient depth at one or more required positions → Indeterminate
- All evaluated positions matching reference alleles with adequate depth → Susceptible

---

# erm41 Truncation Detection

The pipeline evaluates the large deletion associated with truncated `erm41`
using read-depth patterns across the gene.

The canonical truncation contains deletions at approximately:

```text
Δ64–65
Δ159–432
```

Coverage-based detection focuses on the large `Δ159–432` region.

A truncation is inferred when:

```text
median depth in deletion region / median flank depth ≤ 0.10
```

with sufficient coverage in the surrounding region.

The minimum flank coverage threshold is:

```text
20×
```

If surrounding coverage is insufficient, truncation status is reported as
indeterminate rather than assuming the gene is intact.

---

# erm39 Mapping Callability

Because `erm39` homologs may be divergent from the ATCC 35855 reference,
initiation-codon interpretation requires adequate mapping across the gene.

Current requirements are:

```text
≥90% of erm39 bases covered at ≥10×
```

and:

```text
initiation-codon target depth ≥30×
```

If whole-gene mapping does not meet the callability requirement, sequence
variation at the initiation codon is not used as definitive evidence for the
susceptibility-associated genotype.

---

# Output Files

Final aggregate outputs are written to:

```text
<outdir>/summary/
```

The two primary outputs are:

```text
interpretation.tsv
myco_prediction_summary.xlsx
```

`interpretation.tsv` contains the final WGS-based resistance interpretations
and the key genotype evidence used by the interpretation logic.

`myco_prediction_summary.xlsx` provides the underlying resistance evidence in
a convenient review format and can be used to evaluate the evidence supporting
each final interpretation.

Additional TSV files are retained to provide detailed and reproducible evidence
from individual components of the analysis.

---

# Primary Outputs

## 1. `summary/interpretation.tsv`

**Purpose:** Primary machine-readable output containing the final WGS-based
clarithromycin and amikacin resistance interpretations.

Each isolate is represented by a single row containing the final calls and the
key resistance evidence used by the interpretation logic.

| Column | Description |
|---|---|
| Isolate | Isolate identifier |
| clarithromycin_call | Final clarithromycin interpretation |
| clarithromycin_reason | Primary rule or evidence supporting the interpretation |
| amikacin_call | Final amikacin interpretation |
| amikacin_reason | Primary rule or evidence supporting the interpretation |
| `rrl_*` | Observed bases or mixed states at evaluated `rrl` positions |
| `erm41_28` | Observed base at `erm41` position 28 |
| `erm41_19` | Observed base at `erm41` position 19 |
| `rrs_*` | Observed bases or mixed states at evaluated `rrs` positions |
| erm41_truncation | `erm41` truncation assessment |
| erm39_detected | Whether `erm39` was detected by BLAST |
| erm39_1 | `erm39` initiation-codon position 1 genotype |
| erm55_detected | Whether `erm55` was detected by BLAST |
| notes | Additional resistance evidence or interpretation notes |

For site-level genotype columns:

```text
A / C / G / T   Observed nucleotide
REF/ALT         Mixed allele
N               Insufficient evidence / uncallable
NA              Locus not applicable or gene not detected
```

For `erm39_1`, values are interpreted as:

```text
G     Intact GTG initiation codon
C     GTG→CTG susceptibility-associated initiation-codon mutation
G/C   Mixed initiation-codon genotype
N     erm39 detected but initiation-codon position is not callable
NA    erm39 not detected
```

The final interpretation should be reviewed together with the corresponding
evidence in `myco_prediction_summary.xlsx`.

---

## 2. `summary/myco_prediction_summary.xlsx`

**Purpose:** Primary review workbook containing the evidence used to support
the final WGS interpretations.

The workbook aggregates evidence from the major analytical components of the
pipeline. It contains the following worksheets:

- `variants` — site-level sequencing depth and variant evidence from `sites_evidence.tsv`
- `truncation` — coverage metrics used to assess `erm41` truncation
- `blast` — BLAST evidence for `erm41`, `erm39`, and `erm55`
- `coverage` — read depth at each predefined target site for every isolate

The `coverage` worksheet provides a compact per-isolate view of depth at the
evaluated `rrl`, `rrs`, and `erm41` positions and at `erm39` position 1. This
allows uncallable or low-depth sites in the final interpretation to be traced
directly to their observed sequencing depth.

The workbook is intended to provide a convenient format for reviewing the
evidence underlying the calls reported in `interpretation.tsv`.

`interpretation.tsv` remains the primary machine-readable interpretation
output, while the Excel workbook provides the associated evidence for review.

---

# Additional Output Files

## 3. `summary/sites_evidence.tsv`

**Purpose:** Detailed sequencing-depth and variant evidence at predefined
resistance-associated positions.

The table contains one row for every predefined `rrl`, `rrs`, and `erm41`
target site for every isolate, including positions that match the reference
and positions with insufficient depth for interpretation. `erm39` position 1
is also included using the targeted ATCC 35855 mapping assessment. This allows
the observed depth underlying low-coverage or uncallable site-level results to
be reviewed directly.

The table includes evaluated positions from:

- `rrl`
- `rrs`
- `erm41`
- `erm39`

| Column | Description |
|---|---|
| Isolate | Isolate identifier |
| Gene | Resistance-associated gene |
| position | Gene-relative 1-based position |
| Depth | Observed read depth at the target position |
| REF | Reference allele |
| ALT | Alternate allele when variant evidence is present |
| QUAL | Variant quality when applicable |
| DP | Variant-call depth when applicable |
| AD | Allele depths when applicable |
| AF | Alternate allele frequency when applicable |

For reference-matching positions without a variant record, the variant-specific
fields (`ALT`, `QUAL`, `DP`, `AD`, and `AF`) may be blank while `Depth` retains
the observed site coverage.

Allele-frequency interpretation:

```text
AF ≥ 0.90          MUT
0.10 ≤ AF < 0.90  MIXED
AF < 0.10          Reference
Depth < 30×        INDETERMINATE
```

---

## 4. `summary/blast_top_hits.tsv`

**Purpose:** Detailed BLAST evidence for detection of `erm41`, `erm39`, and
`erm55`.

Primary detection threshold:

```text
Percent identity ≥ 90%
AND
Query coverage ≥ 90%
```

| Column | Description |
|---|---|
| Isolate | Isolate identifier |
| GeneGroup | `erm41`, `erm39`, or `erm55` |
| QueryID | Query FASTA sequence identifier |
| QueryLen | Query gene length |
| Subject | Assembly contig containing the hit |
| Pident | Percent nucleotide identity |
| AlnLen | Alignment length |
| QcovPct | Percent query coverage |
| Qstart, Qend | Query alignment coordinates |
| Sstart, Send | Subject alignment coordinates |
| Evalue | BLAST E-value |
| Bitscore | BLAST bit score |
| TotalHits | Total detected hits |
| PreferredHits | Hits satisfying the preferred identity and coverage thresholds |

---

## 5. `summary/erm41_truncation_metrics.tsv`

**Purpose:** Detailed coverage metrics used to assess `erm41` truncation.

| Column | Description |
|---|---|
| Isolate | Isolate identifier |
| erm_len | `erm41` gene length |
| erm_median_depth | Median gene depth |
| left_median_depth | Median left-flank depth |
| right_median_depth | Median right-flank depth |
| flank_median_depth | Combined flank depth |
| erm_del_range | Gene-relative deletion region |
| del_median_depth | Median depth within deletion region |
| del_to_flank_ratio | Deletion-region depth divided by flank depth |
| erm41_callable | Whether coverage is sufficient for truncation assessment |

These metrics provide the underlying evidence for the `erm41` truncation
status incorporated into the final interpretation.

---

# Simulated Validation Dataset

The repository includes an in silico framework for testing expected pipeline
behavior using known genotypes.

The validation dataset includes:

- Unmodified ATCC 19977 reference
- Individual `rrl` resistance-associated mutants
- Individual `rrs` resistance-associated mutants
- `erm41` position 19 and position 28 mutants
- Canonical `erm41` truncation
- Mixed allele datasets
- Low-coverage datasets
- ATCC 35855 `erm39` positive control
- ATCC 35855 `erm39` initiation-codon mutant
- Chromosomal `erm55` control
- Plasmid-associated `erm55` control
- Transposon-associated `erm55` control

## Generate Simulated Dataset

```bash
bash bin/run_atcc_simulation.sh \
    --outdir /path/to/output
```

Help and available options can be displayed with:

```bash
bash bin/run_atcc_simulation.sh --help
```

The simulation workflow creates:

```text
/path/to/output/atcc_dataset/
```

including simulated paired-end reads and a linelist for pipeline execution.

## Run Pipeline on Simulated Dataset

```bash
bash bin/submit_ntm_pipeline.sh \
    /path/to/output/atcc_dataset/linelist.tsv \
    /path/to/output/results \
    --simulated /path/to/output
```

The `--simulated` option directs the pipeline to reads located under:

```text
<DIR>/atcc_dataset/reads/<isolate>/<isolate>_R1.trim.fq.gz
<DIR>/atcc_dataset/reads/<isolate>/<isolate>_R2.trim.fq.gz
```

---

# Validation and Interpretation Scope

The pipeline has been evaluated using both in silico controls with known
expected genotypes and isolates with available phenotypic antimicrobial
susceptibility testing.

In silico validation evaluates:

- correct detection of predefined resistance-associated mutations
- reference and mutant alleles
- mixed allele populations
- low-coverage conditions
- `erm41` genotype and truncation
- `erm39` detection and initiation-codon genotype
- `erm55` detection

WGS interpretations are intentionally restricted to resistance determinants
and variants included within the validated interpretation rules.

A WGS-based susceptible interpretation therefore indicates that a recognized
resistance mechanism was not identified within the current rule set. It does
not exclude resistance caused by uncharacterized mechanisms or variants whose
genotype-phenotype relationships are not sufficiently established for
interpretation.

This is particularly relevant to `erm39`, for which additional sequence
variation may occur but is not interpreted unless the genotype-phenotype
relationship is sufficiently established and incorporated into the pipeline's
predefined rules.

Low sequencing coverage or otherwise insufficient evidence is reported as
indeterminate rather than being interpreted as susceptible.

---

# Repository Structure

```text
uva-ntm_resistance_pipe/
├── bin/          Pipeline and simulation submission scripts
├── config/       Environment configuration
├── workflow/     SLURM workflow jobs
├── scripts/      Analysis and interpretation scripts
├── simulation/   In silico validation scripts
└── references/   Reference genomes and resistance-gene sequences
```

---

# References

## Reference Sequences

- *M. abscessus* ATCC 19977: CU458896.1
- *M. fortuitum* ATCC 35855: CP194182.1
- `erm39`: AY487229.1
- `erm55`: OQ656455.1, OQ656456.1, OQ656457.1

## Software

The pipeline uses:

- SPAdes
- BWA
- SAMtools
- BCFtools
- BLAST+
- Python
