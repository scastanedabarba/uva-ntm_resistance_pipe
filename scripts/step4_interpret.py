#!/usr/bin/env python3
import argparse
import os
import re
from typing import Dict, Tuple, List, Optional

# -----------------------------
# Coordinates (must match step2)
# -----------------------------
RRS_START, RRS_END = 1462398, 1463901
RRL_START, RRL_END = 1464208, 1467319
ERM41_START, ERM41_END = 2345955, 2346476

GENE_STARTS = {"rrs": RRS_START, "rrl": RRL_START, "erm41": ERM41_START}

RRL_SITES = [2269, 2270, 2271, 2281, 2293]
RRS_SITES = [1373, 1375, 1376, 1458]
ERM_SITES = [19, 28]

ALL_SITES = (
    [("rrl", p) for p in RRL_SITES] +
    [("erm41", p) for p in ERM_SITES] +
    [("rrs", p) for p in RRS_SITES]
)

# Truncation rule
TRUNC_RATIO_CUTOFF = 0.10

# erm gene detection notes (as requested)
ERM39_NOTE = "erm39 detected; may confer macrolide resistance in M. fortuitum"
ERM55_NOTE = "erm55 detected; may confer macrolide resistance in M. chelonae"


# -----------------------------
# Small utils
# -----------------------------
def safe_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")

def safe_int(x: str) -> int:
    try:
        return int(x)
    except Exception:
        return -1

def read_tsv(path: str) -> List[List[str]]:
    rows: List[List[str]] = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return rows
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))
    return rows

def read_depth_tsv(depth_tsv: str, ref_contig: str) -> Dict[int, int]:
    """
    samtools depth output: CHROM POS DEPTH
    """
    d: Dict[int, int] = {}
    if not os.path.exists(depth_tsv) or os.path.getsize(depth_tsv) == 0:
        return d
    with open(depth_tsv, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if chrom != ref_contig:
                continue
            pos = safe_int(parts[1])
            dep = safe_int(parts[2])
            if pos >= 1:
                d[pos] = dep if dep >= 0 else 0
    return d

def fasta_load_first_contig(path: str) -> Tuple[str, str]:
    """
    Load the first contig sequence from a FASTA.
    Returns (contig_name, sequence_upper).
    """
    name: Optional[str] = None
    seq_parts: List[str] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is None:
                    name = line[1:].split()[0]
                else:
                    break
            else:
                if name is not None:
                    seq_parts.append(line)
    if name is None:
        raise RuntimeError(f"Could not parse FASTA header from {path}")
    seq = "".join(seq_parts).upper()
    return name, seq

def ref_base(seq: str, ref_pos_1based: int) -> str:
    i = ref_pos_1based - 1
    if i < 0 or i >= len(seq):
        return "N"
    b = seq[i].upper()
    return b if b in ("A", "C", "G", "T") else "N"

def site_ref_pos(gene: str, pos_gene: int) -> int:
    return GENE_STARTS[gene] + (pos_gene - 1)


# -----------------------------
# Parse inputs
# -----------------------------
def parse_sites_evidence(path: str) -> Dict[Tuple[str, str, int], dict]:
    """
    sites_evidence.tsv columns (current):
    Isolate Gene position Depth REF ALT QUAL DP AD AF
    Returns dict keyed by (isolate, gene, pos_gene) -> row dict
    """
    out: Dict[Tuple[str, str, int], dict] = {}
    rows = read_tsv(path)
    if not rows:
        return out
    header = rows[0]
    idx = {h: i for i, h in enumerate(header)}
    for r in rows[1:]:
        if len(r) < len(header):
            r += [""] * (len(header) - len(r))
        iso = r[idx["Isolate"]]
        gene = r[idx["Gene"]]
        pos = safe_int(r[idx["position"]])
        out[(iso, gene, pos)] = {
            "Depth": safe_int(r[idx.get("Depth", -1)]) if idx.get("Depth", -1) >= 0 else -1,
            "REF": r[idx.get("REF", -1)] if idx.get("REF", -1) >= 0 else "",
            "ALT": r[idx.get("ALT", -1)] if idx.get("ALT", -1) >= 0 else "",
            "DP": r[idx.get("DP", -1)] if idx.get("DP", -1) >= 0 else "",
            "AD": r[idx.get("AD", -1)] if idx.get("AD", -1) >= 0 else "",
            "AF": safe_float(r[idx.get("AF", -1)]) if idx.get("AF", -1) >= 0 else float("nan"),
        }
    return out

def parse_trunc_metrics(path: str) -> Dict[str, dict]:
    """
    erm41_truncation_metrics.tsv columns:
    Isolate ... del_to_flank_ratio erm41_callable
    """
    out: Dict[str, dict] = {}
    rows = read_tsv(path)
    if not rows:
        return out
    header = rows[0]
    idx = {h: i for i, h in enumerate(header)}
    for r in rows[1:]:
        if len(r) < len(header):
            r += [""] * (len(header) - len(r))
        iso = r[idx["Isolate"]]
        ratio = safe_float(r[idx.get("del_to_flank_ratio", -1)]) if idx.get("del_to_flank_ratio", -1) >= 0 else float("nan")
        callable_str = r[idx.get("erm41_callable", -1)] if idx.get("erm41_callable", -1) >= 0 else "FALSE"
        callable_flag = True if callable_str.strip().upper() == "TRUE" else False
        out[iso] = {"ratio": ratio, "callable": callable_flag}
    return out

def parse_blast_top_hits(path: str,
                         pident_cutoff: float,
                         qcov_cutoff: float) -> Dict[Tuple[str, str], bool]:
    """
    blast_top_hits.tsv columns include:
    Isolate GeneGroup Pident QcovPct ...
    We'll call detected if Pident>=cutoff and QcovPct>=cutoff.
    Returns dict keyed by (isolate, gene_group) -> detected(True/False)
    """
    detected: Dict[Tuple[str, str], bool] = {}
    rows = read_tsv(path)
    if not rows:
        return detected
    header = rows[0]
    idx = {h: i for i, h in enumerate(header)}

    for r in rows[1:]:
        if len(r) < len(header):
            r += [""] * (len(header) - len(r))

        iso = r[idx["Isolate"]]
        grp = r[idx["GeneGroup"]]

        if grp not in ("erm39", "erm55", "erm41"):
            continue

        pident = safe_float(r[idx.get("Pident", -1)]) if idx.get("Pident", -1) >= 0 else float("nan")
        qcov = safe_float(r[idx.get("QcovPct", -1)]) if idx.get("QcovPct", -1) >= 0 else float("nan")

        ok = (pident == pident and qcov == qcov and pident >= pident_cutoff and qcov >= qcov_cutoff)
        detected[(iso, grp)] = ok

    return detected


# -----------------------------
# Calling primitives
# -----------------------------
def site_state_from_variant(af: float) -> str:
    if af != af:  # nan
        return "WT"
    if af >= 0.90:
        return "MUT"
    if af >= 0.10:
        return "MIXED"
    return "WT"

def call_rrl_rrs_site(iso: str,
                      gene: str,
                      pos_gene: int,
                      dp_min: int,
                      site_map: Dict[Tuple[str, str, int], dict],
                      depths_ref: Dict[int, int]) -> str:
    """
    Return WT/MUT/MIXED/INDETERMINATE for rrl/rrs.
    """
    refpos = site_ref_pos(gene, pos_gene)
    depth = depths_ref.get(refpos, 0)
    if depth < dp_min:
        return "INDETERMINATE"

    key = (iso, gene, pos_gene)
    if key in site_map:
        af = site_map[key].get("AF", float("nan"))
        return site_state_from_variant(af)

    return "WT"

def call_erm41_base(iso: str,
                    pos_gene: int,
                    dp_min: int,
                    site_map: Dict[Tuple[str, str, int], dict],
                    depths_ref: Dict[int, int],
                    ref_seq: str) -> str:
    """
    Return C/T or MIXED or N for erm41 positions (19,28).
    """
    refpos = site_ref_pos("erm41", pos_gene)
    depth = depths_ref.get(refpos, 0)
    if depth < dp_min:
        return "N"

    key = (iso, "erm41", pos_gene)
    if key in site_map:
        af = site_map[key].get("AF", float("nan"))
        st = site_state_from_variant(af)
        if st == "MIXED":
            return "MIXED"
        if st == "MUT":
            alt = (site_map[key].get("ALT", "") or "").strip().upper()
            alt1 = re.split(r"[,\s]+", alt)[0] if alt else ""
            return alt1 if alt1 in ("A", "C", "G", "T") else "N"
        return ref_base(ref_seq, refpos)

    return ref_base(ref_seq, refpos)

def call_erm41_truncation(iso: str, trunc_map: Dict[str, dict]) -> str:
    if iso not in trunc_map:
        return "INDETERMINATE"
    if not trunc_map[iso]["callable"]:
        return "INDETERMINATE"
    ratio = trunc_map[iso]["ratio"]
    if ratio == ratio and ratio <= TRUNC_RATIO_CUTOFF:
        return "TRUNCATED"
    return "NOT_TRUNCATED"


# -----------------------------
# Interpretation
# -----------------------------
def interpret_clarithro(iso: str,
                        rrl_calls: Dict[int, str],
                        erm41_trunc: str,
                        erm41_19: str,
                        erm41_28: str) -> Tuple[str, str]:
    """
    Decision tree:
      - rrl MUT => Resistant
      - rrl MIXED => evaluate erm41 if applicable; else Resistance possible (rrl mixed)
      - rrl WT => evaluate erm41 truncation / bases if applicable; otherwise do not force Indeterminate
    """
    # rrl dominates: any MUT
    for pos in RRL_SITES:
        if rrl_calls.get(pos) == "MUT":
            return ("Resistant", f"Resistant; mutation at rrl position {pos}")

    any_rrl_mixed_pos = next((p for p in RRL_SITES if rrl_calls.get(p) == "MIXED"), None)
    any_rrl_indet = any(rrl_calls.get(p) == "INDETERMINATE" for p in RRL_SITES)

    def erm41_resistance_evidence() -> Optional[str]:
        # Only meaningful if NOT_TRUNCATED and bases callable
        if erm41_trunc != "NOT_TRUNCATED":
            return None
        if erm41_28 == "MIXED":
            return "MIXED_28"
        if erm41_19 == "MIXED":
            return "MIXED_19"
        if erm41_28 == "N" or erm41_19 == "N":
            return "INDET"
        if erm41_28 == "C":
            return None
        if erm41_28 == "T":
            if erm41_19 == "C":
                return "RES_19"
            if erm41_19 == "T":
                return None
        return None

    # If rrl MIXED, still check erm41 (if applicable); otherwise fall back to rrl mixed
    if any_rrl_mixed_pos is not None:
        ev = erm41_resistance_evidence()
        if ev == "RES_19":
            return ("Resistant", "Resistant; erm41 genotype 28=T and 19=C")
        if ev == "MIXED_28":
            return ("Resistance possible", "Resistance possible; mixed genotype at erm41 position 28")
        if ev == "MIXED_19":
            return ("Resistance possible", "Resistance possible; mixed genotype at erm41 position 19")
        return ("Resistance possible", f"Resistance possible; mixed genotype at rrl position {any_rrl_mixed_pos}")

    # If rrl not callable anywhere, we keep indeterminate (unless erm41 provides a decisive call)
    if any_rrl_indet:
        if erm41_trunc == "TRUNCATED":
            return ("Susceptible", "Susceptible; erm41 truncated")
        if erm41_trunc == "NOT_TRUNCATED":
            if erm41_28 == "MIXED":
                return ("Resistance possible", "Resistance possible; mixed genotype at erm41 position 28")
            if erm41_28 == "N":
                return ("Indeterminate", "Indeterminate; insufficient depth at erm41 position 28")
            if erm41_28 == "C":
                return ("Susceptible", "Susceptible; erm41 full-length and position 28 = C")
            if erm41_28 == "T":
                if erm41_19 == "MIXED":
                    return ("Resistance possible", "Resistance possible; mixed genotype at erm41 position 19")
                if erm41_19 == "N":
                    return ("Indeterminate", "Indeterminate; insufficient depth at erm41 position 19")
                if erm41_19 == "C":
                    return ("Resistant", "Resistant; erm41 genotype 28=T and 19=C")
                if erm41_19 == "T":
                    return ("Susceptible", "Susceptible; erm41 full-length with 28=T and 19=T")

        # If erm41 not applicable, don't try to rescue an indeterminate rrl case
        if erm41_trunc == "NOT_APPLICABLE":
            return ("Indeterminate", "Indeterminate; rrl not callable")

        if erm41_trunc == "INDETERMINATE":
            return ("Indeterminate", "Indeterminate; rrl not callable and erm41 not callable for truncation assessment")

        return ("Indeterminate", "Indeterminate; rrl not callable")

    # rrl WT/callable:
    # If erm41 not applicable, do NOT force Indeterminate. (No special reason text requested.)
    if erm41_trunc == "NOT_APPLICABLE":
        return ("Susceptible", "Susceptible")

    if erm41_trunc == "TRUNCATED":
        return ("Susceptible", "Susceptible; erm41 truncated")
    if erm41_trunc == "INDETERMINATE":
        return ("Indeterminate", "Indeterminate; erm41 not callable for truncation assessment")

    # NOT_TRUNCATED
    if erm41_28 == "MIXED":
        return ("Resistance possible", "Resistance possible; mixed genotype at erm41 position 28")
    if erm41_28 == "N":
        return ("Indeterminate", "Indeterminate; insufficient depth at erm41 position 28")
    if erm41_28 == "C":
        return ("Susceptible", "Susceptible; erm41 full-length and position 28 = C")
    if erm41_28 == "T":
        if erm41_19 == "MIXED":
            return ("Resistance possible", "Resistance possible; mixed genotype at erm41 position 19")
        if erm41_19 == "N":
            return ("Indeterminate", "Indeterminate; insufficient depth at erm41 position 19")
        if erm41_19 == "C":
            return ("Resistant", "Resistant; erm41 genotype 28=T and 19=C")
        if erm41_19 == "T":
            return ("Susceptible", "Susceptible; erm41 full-length with 28=T and 19=T")

    return ("Indeterminate", "Indeterminate; unable to resolve erm41 rules")

def interpret_amikacin(iso: str, rrs_calls: Dict[int, str]) -> Tuple[str, str]:
    for pos in RRS_SITES:
        if rrs_calls.get(pos) == "MUT":
            return ("Resistant", f"Resistant; mutation at rrs position {pos}")
    mixed_pos = next((p for p in RRS_SITES if rrs_calls.get(p) == "MIXED"), None)
    if mixed_pos is not None:
        return ("Resistance possible", f"Resistance possible; mixed genotype at rrs position {mixed_pos}")
    if any(rrs_calls.get(p) == "INDETERMINATE" for p in RRS_SITES):
        return ("Indeterminate", "Indeterminate; insufficient depth at one or more rrs target sites")
    return ("Susceptible", "Susceptible; all rrs target sites WT with adequate depth")


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser(description="Step4 interpret: compile per-isolate calls for clarithromycin/amikacin + notes.")
    ap.add_argument("--outdir", required=True, help="Pipeline output root (contains per-isolate subdirs + summary/)")
    ap.add_argument("--ref-fasta", default="", help="Reference FASTA used for mapping (default: <outdir>/scripts/ATCC19977.fasta)")
    ap.add_argument("--ref-contig", default="CU458896.1", help="Reference contig name in targets_depth.tsv (default: CU458896.1)")
    ap.add_argument("--dp-min", type=int, default=30, help="Depth threshold for site callability (default: 30)")
    ap.add_argument("--blast-pident", type=float, default=90.0, help="BLAST pident cutoff for erm detection (default: 90)")
    ap.add_argument("--blast-qcov", type=float, default=90.0, help="BLAST qcov cutoff for erm detection (default: 90)")
    args = ap.parse_args()

    outdir = os.path.abspath(args.outdir)
    summary_dir = os.path.join(outdir, "summary")

    sites_path = os.path.join(summary_dir, "sites_evidence.tsv")
    trunc_path = os.path.join(summary_dir, "erm41_truncation_metrics.tsv")
    blast_path = os.path.join(summary_dir, "blast_top_hits.tsv")

    ref_fasta = args.ref_fasta.strip() or os.path.join(outdir, "scripts", "ATCC19977.fasta")
    if not os.path.exists(ref_fasta):
        raise SystemExit(f"ERROR: Reference FASTA not found: {ref_fasta} (use --ref-fasta to set it)")

    _, ref_seq = fasta_load_first_contig(ref_fasta)

    site_map = parse_sites_evidence(sites_path)
    trunc_map = parse_trunc_metrics(trunc_path)
    blast_detect = parse_blast_top_hits(blast_path, args.blast_pident, args.blast_qcov)

    # Determine isolates to report:
    isolates: List[str] = []
    if trunc_map:
        isolates = sorted(trunc_map.keys())
    else:
        for name in sorted(os.listdir(outdir)):
            p = os.path.join(outdir, name, "variants", "targets_depth.tsv")
            if os.path.exists(p):
                isolates.append(name)

    out_tsv = os.path.join(summary_dir, "interpretation.tsv")

    header = (
        ["Isolate"] +
        ["clarithromycin_call", "clarithromycin_reason",
         "amikacin_call", "amikacin_reason"] +
        [f"rrl_{p}" for p in RRL_SITES] +
        ["erm41_28", "erm41_19"] +
        [f"rrs_{p}" for p in RRS_SITES] +
        ["erm41_truncation"] +
        ["erm39_detected", "erm55_detected"] +
        ["notes"]
    )

    with open(out_tsv, "w") as out:
        out.write("\t".join(header) + "\n")

        for iso in isolates:
            depth_path = os.path.join(outdir, iso, "variants", "targets_depth.tsv")
            depths = read_depth_tsv(depth_path, args.ref_contig)

            # erm detection (BLAST)
            erm41_blast = blast_detect.get((iso, "erm41"), False)
            erm39 = "Y" if blast_detect.get((iso, "erm39"), False) else "N"
            erm55 = "Y" if blast_detect.get((iso, "erm55"), False) else "N"

            # NEW: erm41 applicability gate (coverage primary, BLAST fallback)
            erm41_cov_callable = trunc_map.get(iso, {}).get("callable", False)
            erm41_applicable = bool(erm41_cov_callable or erm41_blast)

            # Call truncation / bases:
            if not erm41_applicable:
                erm41_trunc = "NOT_APPLICABLE"
                erm41_19 = "NA"
                erm41_28 = "NA"
            else:
                erm41_trunc = call_erm41_truncation(iso, trunc_map)
                erm41_19 = call_erm41_base(iso, 19, args.dp_min, site_map, depths, ref_seq)
                erm41_28 = call_erm41_base(iso, 28, args.dp_min, site_map, depths, ref_seq)

            # Call rrl/rrs states
            rrl_calls: Dict[int, str] = {p: call_rrl_rrs_site(iso, "rrl", p, args.dp_min, site_map, depths) for p in RRL_SITES}
            rrs_calls: Dict[int, str] = {p: call_rrl_rrs_site(iso, "rrs", p, args.dp_min, site_map, depths) for p in RRS_SITES}

            # Drug calls
            clar_call, clar_reason = interpret_clarithro(iso, rrl_calls, erm41_trunc, erm41_19, erm41_28)
            amk_call, amk_reason = interpret_amikacin(iso, rrs_calls)

            # Notes + overrides (unchanged)
            notes: List[str] = []
            if erm39 == "Y":
                notes.append(ERM39_NOTE)
                if clar_call != "Resistant":
                    clar_call = "Resistance possible"
                    if not clar_reason.startswith("Resistance possible"):
                        clar_reason = "Resistance possible; erm39 detected"
            if erm55 == "Y":
                notes.append(ERM55_NOTE)
                if clar_call != "Resistant":
                    clar_call = "Resistance possible"
                    if not clar_reason.startswith("Resistance possible"):
                        clar_reason = "Resistance possible; erm55 detected"

            row = []
            row.append(iso)
            row += [clar_call, clar_reason, amk_call, amk_reason]
            row += [rrl_calls[p] for p in RRL_SITES]
            row += [erm41_28, erm41_19]
            row += [rrs_calls[p] for p in RRS_SITES]
            row += [erm41_trunc]
            row += [erm39, erm55]
            row += ["; ".join(notes)]

            out.write("\t".join(row) + "\n")

    print(f"Wrote: {out_tsv}")


if __name__ == "__main__":
    main()

