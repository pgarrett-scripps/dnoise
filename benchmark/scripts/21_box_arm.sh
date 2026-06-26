#!/usr/bin/env bash
# Generate + search the box-centroider arm for Fig 5, reusing the existing
# original/denoised/msms/watershed arms:
#
#   box : dnoise vertical+halo + greedy small-box centroider (--box-centroid)
#
# Transforms all 18 raw .d, then runs one Sage LFQ search (same config/FASTA as the
# main benchmark). Skips work that already exists. dda_5min only.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="$ROOT/data/dda_5min/raw"
BOX="$ROOT/data/dda_5min/box"
FASTA="$ROOT/data/fasta/hybrid.fasta"
CFG="$ROOT/config/sage.json"
DNOISE="$ROOT/../target/release/dnoise"
DCONF="$ROOT/config/dnoise.toml"
SAGE_BIN="$ROOT/tools/sage-0.15.0-beta.1/sage-v0.15.0-beta.1-x86_64-unknown-linux-gnu/sage"
SAGE_BATCH="${SAGE_BATCH:-3}"

mkdir -p "$BOX"

echo "=== [1/2] box-centroid arm: dnoise --box-centroid x18 ==="
for d in "$RAW"/*.d; do
  name="$(basename "$d")"; out="$BOX/$name"
  if [ -d "$out" ]; then echo "skip $name (exists)"; continue; fi
  echo "box: $name"
  "$DNOISE" "$d" "$out" --config "$DCONF" --ms1-polygon --box-centroid >/dev/null
done

search_arm() {  # indir outdir
  local indir="$1" outdir="$2"
  if [ -f "$outdir/lfq.tsv" ]; then echo "skip search (exists): $outdir"; return; fi
  mkdir -p "$outdir"
  local ds=("$indir"/*.d)
  echo "Sage LFQ: ${#ds[@]} runs -> $outdir"
  "$SAGE_BIN" "$CFG" -f "$FASTA" -o "$outdir" --batch-size "$SAGE_BATCH" "${ds[@]}" \
    > "$outdir/sage.log" 2>&1
  echo "  done -> $outdir/lfq.tsv"
}

echo "=== [2/2] Sage search: box ==="
search_arm "$BOX" "$ROOT/results/dda_5min/box"
echo "ALL DONE"
