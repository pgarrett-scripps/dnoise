#!/usr/bin/env bash
# Build + search the two matched intensity-threshold control arms for the UPS DDA
# datasets, for the streak-vs-intensity comparison at MS1 and MS1+MS2 level:
#   intensity       -> MS1-only per-point threshold (T1), MS2 untouched
#                      compare to the streak `denoised` (MS1) arm
#   intensity_msms  -> MS1 (T1) + MS2 (T2) per-point thresholds
#                      compare to the streak `msms` (MS1+MS/MS) arm
# Thresholds are calibrated per dataset by 36_ups_intensity_calib.py (matched to
# the streak filter's MS1/MS2 removal). Both arms apply --ms1-polygon, matching
# the streak arms' gate. Sage searches against the UPS+E.coli DB.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
DNOISE="$(cd .. && pwd)/target/release/dnoise"
UPS_FASTA="$PWD/data/fasta/ups_ecoli.fasta"
SAGE="$PWD/tools/sage-0.15.0-beta.1/sage-v0.15.0-beta.1-x86_64-unknown-linux-gnu/sage"
CFG="$PWD/config/sage.json"

for ds in dda_ups_30spd dda_ups_15spd; do
  raw=$(ls -d data/$ds/raw/*.d)
  echo "######## $ds ########"

  # --- denoise the two intensity arms ---
  iarm="data/$ds/denoised_intensity"
  imarm="data/$ds/denoised_intensity_msms"
  rm -rf "$iarm" "$imarm"
  echo "  [intensity MS1] T1 ..."
  "$DNOISE" "$raw" "$iarm/$(basename "$raw")" \
    --config "config/dnoise.intensity.$ds.toml" --ms1-polygon
  echo "  [intensity MS1+MS2] T1+T2 ..."
  "$DNOISE" "$raw" "$imarm/$(basename "$raw")" \
    --config "config/dnoise.intensity_msms.$ds.toml" --ms1-polygon --denoise-msms

  # --- Sage search both arms (UPS+E.coli DB) ---
  for arm in intensity intensity_msms; do
    out="results/$ds/$arm"
    mkdir -p "$out"
    echo "  [sage $arm] ..."
    "$SAGE" "$CFG" -f "$UPS_FASTA" -o "$out" --batch-size 1 \
      "data/$ds/denoised_$arm"/*.d > "$out/sage.log" 2>&1
  done
done
echo "DONE"
