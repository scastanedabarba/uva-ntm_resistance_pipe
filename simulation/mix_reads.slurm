#!/usr/bin/env bash
#SBATCH --job-name=mix_reads
#SBATCH --output=slurm_logs/mix_reads_%j.out
#SBATCH --error=slurm_logs/mix_reads_%j.err
#SBATCH --account=amr_services
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

: "${WORKDIR:?ERROR: Set WORKDIR env var}"

READSROOT="${WORKDIR%/}/atcc_dataset/reads"

WT_DATASET="${WT_DATASET:-WT}"
MUT_DATASET="${MUT_DATASET:-rrl_2270}"

# Total pairs in the mixed dataset.
# Keep modest by default to avoid huge sampling buffers; override if you want.
MIX_PAIRS="${MIX_PAIRS:-200000}"

SEED="${SEED:-12345}"
FORCE="${FORCE:-0}"

mix_one() {
  local frac="$1"      # 0.05 or 0.50
  local outname="$2"   # rrl_2270_mix05 etc.

  local wt_dir="$READSROOT/$WT_DATASET"
  local mut_dir="$READSROOT/$MUT_DATASET"
  [[ -d "$wt_dir" ]] || { echo "ERROR: WT reads dir missing: $wt_dir" >&2; exit 1; }
  [[ -d "$mut_dir" ]] || { echo "ERROR: MUT reads dir missing: $mut_dir" >&2; exit 1; }

  local outdir="$READSROOT/$outname"
  mkdir -p "$outdir"

  local r1_out="$outdir/${outname}_R1.trim.fq.gz"
  local r2_out="$outdir/${outname}_R2.trim.fq.gz"
  local se_out="$outdir/${outname}_SE.trim.fq.gz"

  if [[ "$FORCE" -eq 0 ]] && [[ -s "$r1_out" && -s "$r2_out" && -s "$se_out" ]]; then
    echo "Skipping (exists): $outname"
    return 0
  fi

  local k_mut k_wt
  k_mut=$(python3 - <<PY
import math
print(int(round($MIX_PAIRS * float("$frac"))))
PY
)
  k_wt=$((MIX_PAIRS - k_mut))
  [[ "$k_mut" -ge 0 ]] || k_mut=0
  [[ "$k_wt" -ge 0 ]] || k_wt=0

  echo "Mixing $outname: total_pairs=$MIX_PAIRS mut_pairs=$k_mut wt_pairs=$k_wt"

  python3 - "$wt_dir" "$mut_dir" "$outdir" "$outname" "$k_wt" "$k_mut" "$SEED" <<'PY'
import gzip, os, random, sys
from pathlib import Path

wt_dir = Path(sys.argv[1])
mut_dir = Path(sys.argv[2])
outdir = Path(sys.argv[3])
outname = sys.argv[4]
k_wt = int(sys.argv[5])
k_mut = int(sys.argv[6])
seed = int(sys.argv[7])

def fq_paths(d: Path, name: str):
    r1 = d / f"{name}_R1.trim.fq.gz"
    r2 = d / f"{name}_R2.trim.fq.gz"
    if not r1.exists() or not r2.exists():
        raise SystemExit(f"Missing inputs: {r1} or {r2}")
    return r1, r2

wt_r1, wt_r2 = fq_paths(wt_dir, wt_dir.name)
mut_r1, mut_r2 = fq_paths(mut_dir, mut_dir.name)

rng = random.Random(seed)

def reservoir_sample_pairs(r1_path: Path, r2_path: Path, k: int):
    # Returns list of (r1_rec, r2_rec) where each rec is 4 lines joined.
    if k == 0:
        return []
    res = []
    n = 0
    with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2:
        while True:
            r1 = [f1.readline() for _ in range(4)]
            r2 = [f2.readline() for _ in range(4)]
            if not r1[0] or not r2[0]:
                break
            n += 1
            rec = ("".join(r1), "".join(r2))
            if len(res) < k:
                res.append(rec)
            else:
                j = rng.randrange(n)
                if j < k:
                    res[j] = rec
    if n < k:
        raise SystemExit(f"Requested k={k} pairs but only found n={n} in {r1_path}")
    return res

wt_pairs = reservoir_sample_pairs(wt_r1, wt_r2, k_wt)
mut_pairs = reservoir_sample_pairs(mut_r1, mut_r2, k_mut)

all_pairs = wt_pairs + mut_pairs
rng.shuffle(all_pairs)

out_r1 = outdir / f"{outname}_R1.trim.fq.gz"
out_r2 = outdir / f"{outname}_R2.trim.fq.gz"
out_se = outdir / f"{outname}_SE.trim.fq.gz"

with gzip.open(out_r1, "wt") as o1, gzip.open(out_r2, "wt") as o2:
    for r1_rec, r2_rec in all_pairs:
        o1.write(r1_rec)
        o2.write(r2_rec)

# Make SE non-empty (copy first FASTQ record from R1 output, consistent with your framework)
with gzip.open(out_r1, "rt") as inp, gzip.open(out_se, "wt") as outp:
    for _ in range(4):
        line = inp.readline()
        if not line:
            break
        outp.write(line)

print(f"Wrote: {out_r1}")
print(f"Wrote: {out_r2}")
print(f"Wrote: {out_se}")
PY
}

mix_one 0.05 "${MUT_DATASET}_mix05"
mix_one 0.50 "${MUT_DATASET}_mix50"

echo "Done mixing."

