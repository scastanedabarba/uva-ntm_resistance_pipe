#!/usr/bin/env bash

# Shared environment setup for UVA NTM resistance pipeline.
#
# Source this file from pipeline/controller/SLURM scripts:
#   source "$REPO_ROOT/config/setup_environment.sh"
#
# This file centralizes all software versions used by the pipeline. It is
# intentionally inert when sourced: scripts call the setup_* function matching
# the tools they need.

AMR_MODULEFILES="/project/amr_services/modulefiles"
NTM_ENV="/project/amr_services/.conda/ntm_resistance_pipe"
NTM_PYTHON="$NTM_ENV/bin/python"

# Pinned project modules
SPADES_MODULE="spades/4.2.0"
BLAST_MODULE="blast+/2.11.0"
BWA_MODULE="bwa/0.7.17"
HTSLIB_MODULE="htslib/1.23"
SAMTOOLS_MODULE="samtools/1.23"
BCFTOOLS_MODULE="bcftools/1.23"

_ntm_die() {
    echo "ERROR: $*" >&2
    return 1 2>/dev/null || exit 1
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || _ntm_die "Required command not found: $1"
}

reset_ntm_modules() {
    command -v module >/dev/null 2>&1 || _ntm_die "Environment Modules/Lmod command is not available"
    [[ -d "$AMR_MODULEFILES" ]] || _ntm_die "AMR Services module path not found: $AMR_MODULEFILES"

    module purge
    module use "$AMR_MODULEFILES"
}

setup_assembly_blast_env() {
    reset_ntm_modules

    # SPAdes loads GCC and Miniforge, which changes the Lmod hierarchy on
    # Rivanna. Re-add the project module path afterward so AMR Services
    # modules retain priority over similarly named system modules.
    module load "$SPADES_MODULE"
    module use "$AMR_MODULEFILES"
    module load "$BLAST_MODULE"

    require_cmd spades.py
    require_cmd makeblastdb
    require_cmd blastn
}

setup_mapping_env() {
    reset_ntm_modules

    module load "$BWA_MODULE"
    module load "$SAMTOOLS_MODULE"
    module load "$BCFTOOLS_MODULE"

    require_cmd bwa
    require_cmd samtools
    require_cmd bcftools
}

setup_mutant_env() {
    reset_ntm_modules

    module load "$SAMTOOLS_MODULE"
    module load "$BCFTOOLS_MODULE"
    # Load HTSlib last so bgzip/tabix resolve from the dedicated HTSlib module.
    module load "$HTSLIB_MODULE"

    require_cmd samtools
    require_cmd bcftools
    require_cmd bgzip
    require_cmd tabix
}

setup_simulation_env() {
    reset_ntm_modules

    module load "$SAMTOOLS_MODULE"

    require_cmd wgsim
}

setup_python_env() {
    [[ -x "$NTM_PYTHON" ]] || _ntm_die "NTM Python environment not found: $NTM_ENV"

    # Prevent ~/.local Python packages from contaminating the fixed environment.
    export PYTHONNOUSERSITE=1

    "$NTM_PYTHON" - <<'PY' >/dev/null 2>&1 || _ntm_die "Required Python packages are not functional in $NTM_ENV"
from Bio import SeqIO
import numpy
import pandas
import matplotlib
import openpyxl
PY
}
