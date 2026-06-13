#!/usr/bin/env bash
# DIA-NN search for a diaPASEF DATASET (dia_5min / dia_15min), 3 arms:
#   original  (RAW)   -> results/<DATASET>/original
#   denoised  (DEN)   -> results/<DATASET>/denoised   (MS1-only dnoise)
#   msms      (MSMS)  -> results/<DATASET>/msms        (MS1 + whole-frame MS2 dnoise)
#
# Engine: DIA-NN (library-free / predicted spectral library), reads Bruker .d
# directly. Set DIANN_BIN to the diann-linux CLI (default: the installed 2.2.0).
#
# SPEED / FAIRNESS: the deep-learning library prediction from the FASTA is the
# expensive step and is acquisition-independent, so we PREDICT IT ONCE from the
# shared hybrid.fasta and REUSE the cached .speclib across every arm and dataset.
# This is what makes the comparison clean -- all arms search the *same* library,
# so the only thing that differs between arms is the (denoised) data. It also
# cuts the predict step from 3x (per arm) to 1x. We deliberately do NOT restrict
# the FASTA to DDA-identified proteins: that would cap DIA's sensitivity at DDA's,
# change the FDR statistics, and could mask the very denoise effect we measure.
#
# Outputs per arm: report.parquet (long) + report.pg_matrix.tsv +
# report.pr_matrix.tsv + report.stats.tsv, read by scripts/12_analyze_dia.py.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/_dataset.sh"
echo "dataset: $DATASET"

DIANN_BIN="${DIANN_BIN:-/home/patrick-garrett/tools/diann-2.2.0/diann-linux}"
THREADS="${DIANN_THREADS:-$(nproc)}"

if [ ! -x "$DIANN_BIN" ]; then
  echo "DIA-NN binary not found / not executable: $DIANN_BIN" >&2
  echo "Set DIANN_BIN=/path/to/diann-linux (see README)." >&2
  exit 2
fi

# DIA-NN loads its bundled shared libraries (libtimsdata.so, libtorch_cpu.so, ...)
# from its own directory; put that on the loader path.
export LD_LIBRARY_PATH="$(dirname "$DIANN_BIN"):${LD_LIBRARY_PATH:-}"
diann() { "$DIANN_BIN" "$@"; }

# DIA-NN needs standard UniProt headers (>sp|ACC|ENTRY). Our shared hybrid.fasta
# carries a SPECIES_ prefix (>HUMAN_sp|ACC|ENTRY) that the Sage path relies on,
# but it breaks DIA-NN's accession/gene parser and collapses the protein list
# (e.g. 31437 -> 6615). Derive a stripped copy for DIA-NN; species are recovered
# downstream in 12_analyze_dia.py by mapping the accession back through the
# original prefixed FASTA (the authoritative tag, so contaminants stay excluded).
DEFAULT_DIANN_FASTA="$ROOT/data/fasta/hybrid_diann.fasta"
DIANN_FASTA="${DIANN_FASTA:-$DEFAULT_DIANN_FASTA}"
# Only auto-derive the default (full proteome) FASTA. A DIANN_FASTA override
# (e.g. the DDA-detected allowlist from scripts/_dda_allowlist.py, a speed win
# that keeps the cross-arm comparison valid) is used as-is.
if [ "$DIANN_FASTA" = "$DEFAULT_DIANN_FASTA" ] \
   && { [ ! -f "$DIANN_FASTA" ] || [ "$FASTA" -nt "$DIANN_FASTA" ]; }; then
  echo "deriving DIA-NN FASTA (standard headers) -> $DIANN_FASTA"
  sed -E 's/^>[A-Za-z]+_(sp|tr)\|/>\1|/' "$FASTA" > "$DIANN_FASTA"
fi
[ -f "$DIANN_FASTA" ] || { echo "DIANN_FASTA not found: $DIANN_FASTA" >&2; exit 1; }
echo "DIA-NN FASTA: $DIANN_FASTA"

# Predicted library is cached PER-FASTA (full vs restricted don't clobber).
LIBDIR="$ROOT/data/fasta/diann_lib/$(basename "${DIANN_FASTA%.fasta}")"
mkdir -p "$LIBDIR"

# In-silico digest + precursor/fragment space for library PREDICTION (only used
# at lib-gen; the search arms just consume the resulting library).
DIGEST=(
  --cut "K*,R*" --missed-cleavages 1
  --min-pep-len 7 --max-pep-len 30
  --min-pr-charge 2 --max-pr-charge 4
  --min-pr-mz 300 --max-pr-mz 1800
  --min-fr-mz 200 --max-fr-mz 1800
  --met-excision --unimod4 --var-mods 1
)

# resolve_speclib: echo the newest *.speclib in LIBDIR, or empty if none.
# `|| true` so a no-match (ls fails under pipefail) doesn't abort the script.
resolve_speclib() { ls -t "$LIBDIR"/*.speclib 2>/dev/null | head -1 || true; }

# 1. Predict the spectral library ONCE from the hybrid FASTA (cached + reused).
SPECLIB="$(resolve_speclib)"
if [ -z "$SPECLIB" ]; then
  echo "=== predicting spectral library from $(basename "$FASTA") (one-time) ==="
  echo "    log: $LIBDIR/gen_lib.log"
  diann \
    --fasta "$DIANN_FASTA" --fasta-search --gen-spec-lib --predictor \
    "${DIGEST[@]}" --threads "$THREADS" \
    --out-lib "$LIBDIR/hybrid.tsv" \
    > "$LIBDIR/gen_lib.log" 2>&1
  SPECLIB="$(resolve_speclib)"
  if [ -z "$SPECLIB" ]; then
    echo "library prediction produced no .speclib -- see $LIBDIR/gen_lib.log:" >&2
    tail -20 "$LIBDIR/gen_lib.log" >&2
    exit 1
  fi
fi
echo "spectral library: $SPECLIB"

# 2. Search each arm against the cached library (no re-prediction).
run_arm() {
  local indir="$1" outdir="$2"
  shopt -s nullglob
  local ds=("$indir"/*.d)
  if [ "${#ds[@]}" -eq 0 ]; then
    echo "skip $(basename "$outdir"): no .d in $indir"
    return 0
  fi
  mkdir -p "$outdir"
  echo "=== DIA-NN [$(basename "$outdir")] : ${#ds[@]} runs -> $outdir ==="
  local args=()
  for d in "${ds[@]}"; do args+=(--f "$d"); done
  diann "${args[@]}" \
    --lib "$SPECLIB" --fasta "$DIANN_FASTA" \
    --reanalyse --matrices --qvalue 0.01 \
    --threads "$THREADS" \
    --out "$outdir/report.parquet" \
    > "$outdir/diann.log" 2>&1
  echo "    done -> $outdir/report.parquet (log: $outdir/diann.log)"
}

[ -d "$RAW" ]  || { echo "no raw .d in $RAW -- run 04_unzip.sh"   >&2; exit 1; }
[ -d "$DEN" ]  || { echo "no MS1-denoised .d in $DEN -- run 04_denoise.sh" >&2; exit 1; }

run_arm "$RAW"  "$RES_ORIGINAL"
run_arm "$DEN"  "$RES_DENOISED"
run_arm "$MSMS" "$RES_MSMS"   # whole-frame MS2 dnoise arm (skipped if absent)
echo "done."
