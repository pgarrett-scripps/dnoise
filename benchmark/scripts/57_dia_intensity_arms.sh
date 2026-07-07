#!/usr/bin/env bash
# Storage-safe build of the two new DIA intensity-threshold control arms
# (intensity, intensity_msms) on the main HYE diaPASEF datasets.
#
# diaPASEF raw data is large (dia_15min raw alone is ~115G for 18 runs) and
# disk is finite, so this NEVER holds more than one new arm's denoised .d
# copies on disk at a time: for each dataset, denoise all 18 runs for ONE arm,
# search it with DIA-NN (full predicted library, fair/denoise-independent --
# see 11_diann.sh), then delete that arm's denoised .d copies once
# report.parquet exists, before moving to the next arm/dataset.
#
# Requires config/dnoise.intensity.<ds>.toml + dnoise.intensity_msms.<ds>.toml
# already calibrated (see 56_recalib_intensity_dia.py).
#
# Usage:
#   bash scripts/57_dia_intensity_arms.sh                      # both gradients
#   DIA_DATASETS="dia_5min" bash scripts/57_dia_intensity_arms.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
DIA_DATASETS="${DIA_DATASETS:-dia_5min dia_15min}"

for ds in $DIA_DATASETS; do
  echo "######## $ds ########"
  DATASET="$ds"; source scripts/_dataset.sh

  for arm in intensity intensity_msms; do
    outdir_var=$([ "$arm" = "intensity" ] && echo INT || echo INTMSMS)
    resdir_var=$([ "$arm" = "intensity" ] && echo RES_INTENSITY || echo RES_INTMSMS)
    outdir="${!outdir_var}"
    resdir="${!resdir_var}"

    if [ -f "$resdir/report.parquet" ]; then
      echo "=== $ds/$arm: report.parquet already exists -- skipping build+search ==="
      continue
    fi

    echo "=== $ds/$arm: denoise (18 runs) ==="
    DATASET="$ds" DENOISE_ARMS="$arm" bash scripts/04_denoise.sh

    echo "=== $ds/$arm: DIA-NN full-lib search ==="
    DATASET="$ds" DIANN_FULL_LIB=1 bash scripts/11_diann.sh

    if [ -f "$resdir/report.parquet" ]; then
      echo "=== $ds/$arm: search complete -- reclaiming disk ($outdir) ==="
      rm -rf "$outdir"
    else
      echo "WARNING: $ds/$arm search did not produce $resdir/report.parquet -- NOT deleting $outdir (inspect $resdir/diann.log)" >&2
    fi
  done
done
echo "ALL DONE"
