#!/usr/bin/env bash
# Unattended end-to-end run for one DATASET:
#   download -> unzip -> denoise -> search -> analyze -> data-reduction.
# Each step is resumable, so re-running continues where it left off.
# Select the dataset with DATASET=dda_15min bash scripts/run_all.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export DATASET="${DATASET:-dda_5min}"
echo "[run_all] $(date) dataset=$DATASET"

echo "[run_all] $(date) === download ==="
bash scripts/03_download.sh

echo "[run_all] $(date) === unzip ==="
bash scripts/04_unzip.sh

echo "[run_all] $(date) === denoise (MS1 + MS/MS arms) ==="
bash scripts/04_denoise.sh

echo "[run_all] $(date) === search (3 arms) ==="
SAGE_BATCH="${SAGE_BATCH:-2}" bash scripts/05_search.sh

echo "[run_all] $(date) === analyze ==="
uv run scripts/06_analyze.py

echo "[run_all] $(date) === data reduction ==="
uv run scripts/07_data_reduction.py

echo "[run_all] $(date) === DONE ==="
