#!/usr/bin/env bash
# Generate + search two extra dda_5min arms for the centroider (Fig 5) and Bruker
# (Fig 6) comparisons, reusing the existing original/denoised/msms/watershed arms:
#
#   box     : dnoise vertical+halo + greedy small-box centroider (--box-centroid)
#   bruker  : Bruker tdf2tdf Minesweeper, --ms1-min-frequency 30 (matched to
#             dnoise's ~32% MS1 retention)
#
# Each arm: transform all 18 raw .d, then one Sage LFQ search (same config/FASTA
# as the main benchmark). Skips work that already exists. dda_5min only.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW="$ROOT/data/dda_5min/raw"
BOX="$ROOT/data/dda_5min/box"
BRUKER="$ROOT/data/dda_5min/bruker"
FASTA="$ROOT/data/fasta/hybrid.fasta"
CFG="$ROOT/config/sage.json"
DNOISE="$ROOT/../target/release/dnoise"
DCONF="$ROOT/config/dnoise.toml"
BRUKER_BIN="${BRUKER_BIN:-/home/patrick-garrett/Downloads/tdf2tdf/tdf-to-tdf_transform}"
BRUKER_DIR="$(dirname "$BRUKER_BIN")"
SAGE_BIN="$ROOT/tools/sage-0.15.0-beta.1/sage-v0.15.0-beta.1-x86_64-unknown-linux-gnu/sage"
SAGE_BATCH="${SAGE_BATCH:-3}"

mkdir -p "$BOX" "$BRUKER"

echo "=== [1/4] box-centroid arm: dnoise --box-centroid x18 ==="
for d in "$RAW"/*.d; do
  name="$(basename "$d")"; out="$BOX/$name"
  if [ -d "$out" ]; then echo "skip $name (exists)"; continue; fi
  echo "box: $name"
  "$DNOISE" "$d" "$out" --config "$DCONF" --box-centroid >/dev/null
done

echo "=== [2/4] bruker arm: tdf2tdf Minesweeper min-freq 30 x18 ==="
for d in "$RAW"/*.d; do
  name="$(basename "$d")"; out="$BRUKER/$name"
  if [ -d "$out" ]; then echo "skip $name (exists)"; continue; fi
  echo "bruker: $name"
  # tdf2tdf writes the result NESTED as <root>/<name>; use a per-file temp root.
  tmp="$BRUKER/.tmp_$name"; rm -rf "$tmp"; mkdir -p "$tmp"
  LD_LIBRARY_PATH="$BRUKER_DIR" "$BRUKER_BIN" -i "$d" -o "$tmp" -c -d \
    --ms1-min-frequency 30 >/dev/null 2>&1
  if [ -d "$tmp/$name" ]; then mv "$tmp/$name" "$out"; else
    echo "ERROR: bruker produced no output for $name" >&2; ls "$tmp" >&2; exit 1
  fi
  rm -rf "$tmp"
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

echo "=== [3/4] Sage search: box ==="
search_arm "$BOX" "$ROOT/results/dda_5min/box"
echo "=== [4/4] Sage search: bruker ==="
search_arm "$BRUKER" "$ROOT/results/dda_5min/bruker"
echo "ALL DONE"
