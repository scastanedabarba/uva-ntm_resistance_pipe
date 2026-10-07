#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# sync_ntm_results.sh
#
# Sync completed UVA NTM resistance pipeline results to:
#
#   /project/amr_services/ntm_resistance/<run>/<isolate>/
#
# The linelist must contain:
#   isolate,run
#
# Files ending in .bam or .bam.bai are excluded.
#
# Existing isolate directories are NEVER overwritten or merged.
# Each isolate is first copied to a temporary directory and only
# renamed to the final isolate name after rsync succeeds.
#
# USAGE:
#   bash sync_ntm_results.sh <linelist.csv|tsv> <workdir>
#
# EXAMPLE:
#   bash sync_ntm_results.sh \
#       /scratch/sgj4qr/erm_pipetest/linelist.csv \
#       /scratch/sgj4qr/erm_pipetest
# ============================================================

usage() {
    cat <<'EOF'
USAGE:
    bash sync_ntm_results.sh <linelist.csv|tsv> <workdir>

ARGUMENTS:
    linelist    CSV or TSV containing isolate and run columns
    workdir     Pipeline output directory containing isolate directories

DESTINATION:
    /project/amr_services/ntm_resistance/<run>/<isolate>/

NOTES:
    - Existing isolate directories are skipped.
    - .bam and .bam.bai files are excluded.
    - Missing source isolate directories are reported and skipped.
    - Transfers use a temporary directory to prevent incomplete
      results from appearing as completed isolate directories.
EOF
}

# ------------------------------------------------------------
# Arguments
# ------------------------------------------------------------

if [[ $# -ne 2 ]]; then
    usage >&2
    exit 1
fi

LINELIST="$1"
WORKDIR="$2"

DEST_ROOT="/project/amr_services/ntm_resistance"

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------

[[ -f "$LINELIST" ]] || {
    echo "ERROR: Linelist not found: $LINELIST" >&2
    exit 1
}

[[ -d "$WORKDIR" ]] || {
    echo "ERROR: Workdir not found: $WORKDIR" >&2
    exit 1
}

command -v rsync >/dev/null 2>&1 || {
    echo "ERROR: rsync is not available." >&2
    exit 1
}

# Normalize paths
LINELIST="$(cd "$(dirname "$LINELIST")" && pwd)/$(basename "$LINELIST")"
WORKDIR="$(cd "$WORKDIR" && pwd)"

mkdir -p "$DEST_ROOT"

# ------------------------------------------------------------
# Parse linelist
#
# Accepts CSV or TSV.
# Header is optional.
# Outputs:
#   isolate<TAB>run
# ------------------------------------------------------------

parse_linelist() {
    awk -v FS='[,\t]' '
        BEGIN {
            OFS="\t"
        }

        NR == 1 {
            if (tolower($1) ~ /isolate/ || tolower($2) ~ /run/) {
                next
            }
        }

        NF >= 2 && $1 != "" && $2 != "" {
            isolate=$1
            run=$2

            gsub(/^[[:space:]]+|[[:space:]]+$/, "", isolate)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", run)

            print isolate, run
        }
    ' "$1"
}

# ------------------------------------------------------------
# Counters
# ------------------------------------------------------------

COPIED=0
SKIPPED=0
MISSING=0
FAILED=0

echo "============================================================"
echo "UVA NTM Resistance Results Sync"
echo "============================================================"
echo "Linelist:    $LINELIST"
echo "Workdir:     $WORKDIR"
echo "Destination: $DEST_ROOT"
echo

# ------------------------------------------------------------
# Sync isolates
# ------------------------------------------------------------

while IFS=$'\t' read -r ISOLATE RUN; do

    [[ -z "$ISOLATE" || -z "$RUN" ]] && continue

    SOURCE="$WORKDIR/$ISOLATE"
    RUN_DEST="$DEST_ROOT/$RUN"
    FINAL_DEST="$RUN_DEST/$ISOLATE"
    TEMP_DEST="$RUN_DEST/.${ISOLATE}.syncing"

    echo "------------------------------------------------------------"
    echo "Isolate: $ISOLATE"
    echo "Run:     $RUN"

    # --------------------------------------------------------
    # Check source
    # --------------------------------------------------------

    if [[ ! -d "$SOURCE" ]]; then
        echo "STATUS: MISSING"
        echo "Source directory not found:"
        echo "  $SOURCE"

        ((MISSING+=1))
        continue
    fi

    # --------------------------------------------------------
    # Never overwrite an existing isolate
    # --------------------------------------------------------

    if [[ -e "$FINAL_DEST" ]]; then
        echo "STATUS: SKIPPED"
        echo "Destination already exists:"
        echo "  $FINAL_DEST"

        ((SKIPPED+=1))
        continue
    fi

    # --------------------------------------------------------
    # Create run directory
    # --------------------------------------------------------

    mkdir -p "$RUN_DEST"

    # --------------------------------------------------------
    # Protect against a stale temporary directory
    # --------------------------------------------------------

    if [[ -e "$TEMP_DEST" ]]; then
        echo "WARNING: Removing stale temporary sync directory:"
        echo "  $TEMP_DEST"
        rm -rf "$TEMP_DEST"
    fi

    mkdir -p "$TEMP_DEST"

    # --------------------------------------------------------
    # Copy isolate
    #
    # Exclude BAM files to avoid unnecessary project storage.
    # The trailing "/" on SOURCE copies the CONTENTS of the
    # isolate directory into TEMP_DEST.
    # --------------------------------------------------------

    echo "Copying:"
    echo "  $SOURCE"
    echo "    ->"
    echo "  $FINAL_DEST"

    if rsync -a \
        --no-owner \
        --no-group \
        --exclude='*.bam' \
        --exclude='*.bam.bai' \
        "$SOURCE/" \
        "$TEMP_DEST/"; then

        # Atomically move completed transfer to final name.
        mv "$TEMP_DEST" "$FINAL_DEST"

        echo "STATUS: COPIED"
        ((COPIED+=1))

    else
        echo "STATUS: FAILED"
        echo "rsync failed for $ISOLATE." >&2

        # Do not leave a partial transfer behind.
        rm -rf "$TEMP_DEST"

        ((FAILED+=1))
    fi

done < <(parse_linelist "$LINELIST")

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Sync Summary"
echo "============================================================"
echo "Copied:  $COPIED"
echo "Skipped: $SKIPPED"
echo "Missing: $MISSING"
echo "Failed:  $FAILED"
echo "============================================================"

if [[ "$FAILED" -gt 0 ]]; then
    exit 1
fi

exit 0
