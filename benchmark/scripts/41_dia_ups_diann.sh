#!/usr/bin/env bash
# DIA-NN streak-vs-intensity benchmark on the UPS2 / timsTOF Pro 2 diaPASEF data
# (the diaPASEF analog of 37_ups_intensity_arms.sh, which did the Sage/DDA arms).
#
# Arms per dataset (5), each searched INDEPENDENTLY against one predicted library:
#   original         raw .d                                   -> results/<ds>/original
#   denoised         streak MS1 (vertical+halo, --dia-ms1-window)
#   intensity        matched per-point MS1 threshold (T1, --dia-ms1-window)
#   msms             streak MS1 + streak MS/MS (--dia-ms1-window --dia-window)
#   intensity_msms   matched per-point MS1 (T1) + MS2 (T2)
# The streak `denoised`/`msms` arms already exist (built by 04_denoise.sh); this
# script builds the two `intensity*` arms (thresholds from 40_..calib.py) and
# runs DIA-NN on all five.
#
# LIBRARY: we predict ONE spectral library from the small UPS+E.coli FASTA and
# search every arm against it independently (DIA-NN --reanalyse). This is the
# diaPASEF analog of the DDA comparison (each Sage arm searched the same FASTA on
# its own) and avoids the raw-library-reuse circularity reviewers flagged: a
# denoised arm can gain or lose IDs on its own merits, not be capped at the raw-
# detectable set. The FASTA is ~6.5x smaller than the HYE proteome, so the one-
# time prediction is cheap relative to the full HYE library.
#
# Idempotent: each arm is skipped if its report.parquet already exists
# (DIANN_FORCE=1 to override); the predicted lib is cached and reused.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
ROOT="$PWD"
DNOISE="$(cd .. && pwd)/target/release/dnoise"
UPS_FASTA="$ROOT/data/fasta/ups_ecoli.fasta"
DIANN_BIN="${DIANN_BIN:-/home/patrick-garrett/tools/diann-2.2.0/diann-linux}"
THREADS="${DIANN_THREADS:-$(nproc)}"
DATASETS="${DIA_UPS_DATASETS:-dia_ups_30spd dia_ups_15spd}"

[ -x "$DIANN_BIN" ] || { echo "DIA-NN binary missing/not executable: $DIANN_BIN" >&2; exit 2; }
[ -x "$DNOISE" ] || { echo "dnoise binary missing: $DNOISE (cargo build --release)" >&2; exit 2; }
export LD_LIBRARY_PATH="$(dirname "$DIANN_BIN"):${LD_LIBRARY_PATH:-}"
diann() { "$DIANN_BIN" "$@"; }

# DIA-NN needs standard >sp|ACC|ENTRY headers; ups_ecoli.fasta carries a SPECIES_
# prefix (>HUMAN_sp|...). Derive a stripped copy (species recovered downstream by
# accession). Accession set is unchanged, so contaminants/UPS/E.coli all survive.
DIANN_FASTA="$ROOT/data/fasta/ups_ecoli_diann.fasta"
if [ ! -f "$DIANN_FASTA" ] || [ "$UPS_FASTA" -nt "$DIANN_FASTA" ]; then
  echo "deriving stripped DIA-NN FASTA -> $(basename "$DIANN_FASTA")"
  sed -E 's/^>[A-Za-z]+_(sp|tr)\|/>\1|/' "$UPS_FASTA" > "$DIANN_FASTA"
fi

# --- predict the spectral library once (cached) ---
LIBDIR="$ROOT/data/fasta/diann_lib/ups_ecoli_diann"
mkdir -p "$LIBDIR"
DIGEST=(
  --cut "K*,R*" --missed-cleavages 1
  --min-pep-len 7 --max-pep-len 30
  --min-pr-charge 2 --max-pr-charge 4
  --min-pr-mz 300 --max-pr-mz 1800
  --min-fr-mz 200 --max-fr-mz 1800
  --met-excision --unimod4 --var-mods 1
)
resolve_speclib() { ls -t "$LIBDIR"/*.speclib 2>/dev/null | head -1 || true; }
SPECLIB="$(resolve_speclib)"
if [ -z "$SPECLIB" ]; then
  echo "=== predicting spectral library from ups_ecoli (one-time) -> $LIBDIR/gen_lib.log ==="
  diann --fasta "$DIANN_FASTA" --fasta-search --gen-spec-lib --predictor \
    "${DIGEST[@]}" --threads "$THREADS" --out-lib "$LIBDIR/ups_ecoli.tsv" \
    > "$LIBDIR/gen_lib.log" 2>&1
  SPECLIB="$(resolve_speclib)"
  [ -n "$SPECLIB" ] || { echo "library prediction produced no .speclib:" >&2; tail -20 "$LIBDIR/gen_lib.log" >&2; exit 1; }
fi
echo "spectral library: $SPECLIB"

run_arm() {  # <indir> <outdir>
  local indir="$1" outdir="$2"
  shopt -s nullglob
  local ds=("$indir"/*.d)
  [ "${#ds[@]}" -gt 0 ] || { echo "skip $(basename "$outdir"): no .d in $indir"; return 0; }
  mkdir -p "$outdir"
  if [ "${DIANN_FORCE:-0}" != "1" ] && [ -f "$outdir/report.parquet" ]; then
    echo "skip (exists): $(basename "$outdir")"; return 0
  fi
  echo "=== DIA-NN [$(basename "$(dirname "$outdir")")/$(basename "$outdir")] : ${#ds[@]} run(s) ==="
  local args=(); for d in "${ds[@]}"; do args+=(--f "$d"); done
  diann "${args[@]}" --lib "$SPECLIB" --fasta "$DIANN_FASTA" \
    --reanalyse --matrices --qvalue 0.01 --threads "$THREADS" \
    --out "$outdir/report.parquet" > "$outdir/diann.log" 2>&1
  echo "    done -> $outdir/report.parquet"
}

for ds in $DATASETS; do
  echo "######## $ds ########"
  raw=$(ls -d data/$ds/raw/*.d)
  rawname=$(basename "$raw")

  # --- build the two matched intensity arms (DIA gates) ---
  icfg="config/dnoise.intensity.$ds.toml"
  imcfg="config/dnoise.intensity_msms.$ds.toml"
  [ -f "$icfg" ] && [ -f "$imcfg" ] || { echo "missing calibrated config(s); run 40_dia_ups_intensity_calib.py" >&2; exit 1; }
  iarm="data/$ds/denoised_intensity"
  imarm="data/$ds/denoised_intensity_msms"
  if [ ! -d "$iarm/$rawname" ] || [ "${DENOISE_FORCE:-0}" = "1" ]; then
    echo "  [intensity MS1] T1 (--dia-ms1-window) ..."
    rm -rf "$iarm"; "$DNOISE" "$raw" "$iarm/$rawname" --config "$icfg" --dia-ms1-window
  else echo "  [intensity MS1] up to date"; fi
  if [ ! -d "$imarm/$rawname" ] || [ "${DENOISE_FORCE:-0}" = "1" ]; then
    echo "  [intensity MS1+MS2] T1+T2 (--dia-ms1-window --dia-window) ..."
    rm -rf "$imarm"; "$DNOISE" "$raw" "$imarm/$rawname" --config "$imcfg" --dia-ms1-window --dia-window --denoise-msms
  else echo "  [intensity MS1+MS2] up to date"; fi

  # --- DIA-NN on all 5 arms ---
  run_arm "data/$ds/raw"                   "results/$ds/original"
  run_arm "data/$ds/denoised"              "results/$ds/denoised"
  run_arm "data/$ds/denoised_intensity"    "results/$ds/intensity"
  run_arm "data/$ds/denoised_msms"         "results/$ds/msms"
  run_arm "data/$ds/denoised_intensity_msms" "results/$ds/intensity_msms"
done
echo "ALL DONE"
