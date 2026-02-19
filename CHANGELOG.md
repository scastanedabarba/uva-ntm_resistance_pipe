UVA NTM Resistance Pipeline CHANGELOG

## [0.2.0] - 2026-02-19

Major Improvements:

-   Refactored pipeline directory structure for clarity and
    reproducibility:

    -   bin/
    -   workflow/
    -   scripts/
    -   simulation/
    -   references/

-   Standardized path handling using REPO_ROOT for consistent script
    referencing.

-   Simplified controller interface to use a single required outdir
    argument.

-   Removed unnecessary isolate argument from summary step.

-   Added –simulated flag to clearly distinguish simulated dataset
    execution from clinical runs.

-   Improved SLURM job structure and dependency handling.

-   Fixed outdir parsing issues that caused summary step directory
    nesting errors.

-   Stabilized SPAdes execution by explicitly setting memory usage.

-   Expanded README with:

    -   Detailed clinical dataset instructions
    -   Clear coordinate definitions (gene-relative 1-based positions)
    -   Explicit site-of-interest listing
    -   Interpretation logic documentation
    -   Reference genome accessions and citations
    -   Simulation workflow documentation

-   Improved output documentation with column-level explanations.

-   Added Excel summary as convenience output (TSV remains primary
    reproducible output).


## [0.1.0] - Initial Release

Features:

-   SPAdes assembly of isolates
-   BLAST detection of erm41, erm39, erm55
-   Mapping to ATCC 19977 reference (CU458896.1)
-   Variant calling at resistance-associated loci
-   erm41 truncation detection (coverage-based)
-   Structured interpretation for:
    -   Clarithromycin
    -   Amikacin
-   Generated reproducible TSV outputs:
    -   blast_top_hits.tsv
    -   sites_evidence.tsv
    -   erm41_truncation_metrics.tsv
    -   interpretation.tsv
-   Designed for SLURM-based execution on UVA Rivanna HPC.

