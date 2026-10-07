#!/usr/bin/env bash
# Shared provenance helpers for the UVA NTM resistance pipeline.
# This file is sourced by workflow jobs; it does not alter analytical behavior.

provenance_init() {
    local iso_dir="$1"
    mkdir -p "$iso_dir/.provenance"
}

provenance_log() {
    local log_file="$1"
    local command_text="$2"
    printf '[%s]\n%s\n\n' "$(date '+%Y-%m-%d %H:%M:%S %z')" "$command_text" >> "$log_file"
}

provenance_section() {
    local log_file="$1"
    local title="$2"
    {
        printf '============================================================\n'
        printf '%s\n' "$title"
        printf '============================================================\n\n'
    } >> "$log_file"
}

provenance_version_line() {
    local out_file="$1"
    local label="$2"
    shift 2
    local value
    value="$("$@" 2>&1 | head -n 1 || true)"
    [[ -n "$value" ]] || value="UNAVAILABLE"
    printf '%s: %s\n' "$label" "$value" >> "$out_file"
}

