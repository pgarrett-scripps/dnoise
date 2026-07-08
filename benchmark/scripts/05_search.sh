#!/usr/bin/env bash
# Run Sage on the three arms of a dataset with one shared config + FASTA:
#   original -> results/<DATASET>/original   (data/<DATASET>/raw)
#   denoised -> results/<DATASET>/denoised   (data/<DATASET>/denoised, MS1)
#   msms     -> results/<DATASET>/msms        (data/<DATASET>/denoised_msms)
#   intensity / intensity_msms -> matched intensity-threshold control arms
# Each optional arm is searched only if its .d folders exist.
#
# Uses the Sage 0.15.0-beta.1 prebuilt, which does timsTOF ion-mobility LFQ
# directly on .d (verified: populated lfq.tsv). The 0.14.6 stable that is usually
# on PATH writes an EMPTY lfq.tsv (no mobility LFQ), and the 0.15.0-beta.2
# prebuilt panics resolving the .d path (timsrust 0.4.2 `with_path`), so we pin
# beta.1. Override with SAGE_BIN=/path/to/sage if you have your own build.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/_dataset.sh"
CFG="$ROOT/config/sage.json"
BATCH="${SAGE_BATCH:-2}"   # files searched in parallel; lower if RAM-limited
echo "dataset: $DATASET"

SAGE_VER="0.15.0-beta.1"
SAGE_DIR="$ROOT/tools/sage-$SAGE_VER"
SAGE_BIN="${SAGE_BIN:-$SAGE_DIR/sage-v$SAGE_VER-x86_64-unknown-linux-gnu/sage}"

# Fetch the Linux x86_64 prebuilt on first use (set SAGE_BIN on other platforms).
if [ ! -x "$SAGE_BIN" ]; then
  echo "fetching Sage $SAGE_VER prebuilt -> $SAGE_DIR"
  mkdir -p "$SAGE_DIR"
  url="https://github.com/lazear/sage/releases/download/v$SAGE_VER/sage-v$SAGE_VER-x86_64-unknown-linux-gnu.tar.gz"
  curl -fsSL "$url" | tar xz -C "$SAGE_DIR"
  chmod +x "$SAGE_BIN"
fi
echo "using $("$SAGE_BIN" --version)  ($SAGE_BIN)"

run_arm () {
  local out="$1"; shift
  mkdir -p "$out"
  # Resumable: skip a completed arm so a re-run after interruption continues
  # instead of redoing it. Set SAGE_FORCE=1 to re-search regardless.
  if [ "${SAGE_FORCE:-0}" != "1" ] && [ -f "$out/results.sage.tsv" ] && [ -f "$out/lfq.tsv" ]; then
    echo "skip (exists): $(basename "$out")"
    return 0
  fi
  echo "=== Sage [$(basename "$out")] -> $out ==="
  "$SAGE_BIN" "$CFG" -f "$FASTA" -o "$out" --batch-size "$BATCH" "$@"
}

shopt -s nullglob
raw=( "$RAW"/*.d )
den=( "$DEN"/*.d )
msms=( "$MSMS"/*.d )
wshed=( "$WSHED"/*.d )
intensity=( "$INT"/*.d )
intmsms=( "$INTMSMS"/*.d )
[ ${#raw[@]} -gt 0 ] || { echo "no raw .d in $RAW — run 03/04_unzip"; exit 1; }

run_arm "$RES_ORIGINAL" "${raw[@]}"
# The streak MS1 arm's .d copies are routinely deleted after its own search
# completes (storage discipline for later arm additions), so a missing $DEN is
# only a hard error if there's no prior search to fall back on -- run_arm's own
# "skip (exists)" check handles the already-searched case even with den empty.
if [ ${#den[@]} -gt 0 ]; then
  run_arm "$RES_DENOISED" "${den[@]}"
elif [ -f "$RES_DENOISED/results.sage.tsv" ] && [ -f "$RES_DENOISED/lfq.tsv" ]; then
  echo "skip (exists, .d already cleaned up): denoised"
else
  echo "no MS1-denoised .d in $DEN and no prior results in $RES_DENOISED — run 04_denoise.sh" >&2
  exit 1
fi
if [ ${#msms[@]} -gt 0 ]; then
  run_arm "$RES_MSMS" "${msms[@]}"
else
  echo "no MS/MS-denoised .d in $MSMS — skipping msms arm"
fi
if [ ${#wshed[@]} -gt 0 ]; then
  run_arm "$RES_WSHED" "${wshed[@]}"
else
  echo "no watershed-centroided .d in $WSHED — skipping watershed arm"
fi
if [ ${#intensity[@]} -gt 0 ]; then
  run_arm "$RES_INTENSITY" "${intensity[@]}"
else
  echo "no intensity-threshold .d in $INT — skipping intensity arm"
fi
if [ ${#intmsms[@]} -gt 0 ]; then
  run_arm "$RES_INTMSMS" "${intmsms[@]}"
else
  echo "no intensity_msms .d in $INTMSMS — skipping intensity_msms arm"
fi
echo "done."
