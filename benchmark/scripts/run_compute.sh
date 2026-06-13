#!/usr/bin/env bash
# Compute-only pipeline for one DATASET (denoise -> search -> analyze).
# Select the dataset with DATASET=dda_15min bash scripts/run_compute.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
export DATASET="${DATASET:-dda_5min}"
echo "[compute] $(date) dataset=$DATASET"
echo "[compute] $(date) === denoise (MS1 + MS/MS arms) ==="
bash scripts/04_denoise.sh
echo "[compute] $(date) === search (3 arms) ==="
SAGE_BATCH="${SAGE_BATCH:-2}" bash scripts/05_search.sh
echo "[compute] $(date) === analyze ==="
uv run scripts/06_analyze.py
echo "[compute] $(date) === data reduction ==="
uv run scripts/07_data_reduction.py
echo "[compute] $(date) === DONE ==="
