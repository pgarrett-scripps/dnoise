#!/usr/bin/env bash
# Single-sample denoiser intensity comparison.
#
# Searches ONE acquisition (default Condition_A_REP1), denoised seven different
# ways, in a SINGLE Sage LFQ run -- so MS1 features are matched across methods
# (MBR) and the resulting lfq.tsv has one intensity column per method, all for
# the SAME underlying peptides. That isolates how each denoiser shifts the
# quantified intensity of identical peptides.
#
# Arms (all derived from $FILE):
#   original, dnoise_ms1, dnoise_msms,
#   bruker_minesweeper, bruker_minesweeper_strong, bruker_eh, bruker_background
#
# All seven .d share a basename and Sage keys LFQ columns by filename, so we
# stage uniquely-named symlinks (stage/<label>.d). Intensities are deliberately
# NOT normalized downstream: same sample, so any difference IS the denoising
# effect. Analyze with scripts/19_single_compare_analyze.py.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FILE="${FILE:-LFQ_Ultra2_PASEF_5min_50ng_Condition_A_REP1.d}"
FASTA="$ROOT/data/fasta/hybrid.fasta"
CFG="$ROOT/config/sage.json"
OUT="$ROOT/results/single_compare"
STAGE="$OUT/stage"

# Sage 0.15.0-beta.1 (mobility-aware LFQ); same discovery as 05_search.sh.
SAGE_VER="0.15.0-beta.1"
SAGE_DIR="$ROOT/tools/sage-$SAGE_VER"
SAGE_BIN="${SAGE_BIN:-$SAGE_DIR/sage-v$SAGE_VER-x86_64-unknown-linux-gnu/sage}"
if [ ! -x "$SAGE_BIN" ]; then
  echo "fetching Sage $SAGE_VER prebuilt -> $SAGE_DIR"
  mkdir -p "$SAGE_DIR"
  url="https://github.com/lazear/sage/releases/download/v$SAGE_VER/sage-v$SAGE_VER-x86_64-unknown-linux-gnu.tar.gz"
  curl -fsSL "$url" | tar xz -C "$SAGE_DIR"
  chmod +x "$SAGE_BIN"
fi
echo "using $("$SAGE_BIN" --version)"
echo "file:  $FILE"

D5="$ROOT/data/dda_5min"
BR="$D5/bruker_demo"
# label -> source .d  (Bruker outputs are nested one level under the mode dir)
declare -A ARMS=(
  [original]="$D5/raw/$FILE"
  [dnoise_ms1]="$D5/denoised/$FILE"
  [dnoise_msms]="$D5/denoised_msms/$FILE"
  [bruker_minesweeper]="$BR/minesweeper/$FILE"
  [bruker_minesweeper_mf30]="$BR/ms_mf30/$FILE"
  [bruker_minesweeper_strong]="$BR/minesweeper_strong/$FILE"
  [bruker_eh]="$BR/eh/$FILE"
  [bruker_background]="$BR/background/$FILE"
)
ORDER=(original dnoise_ms1 dnoise_msms bruker_minesweeper bruker_minesweeper_mf30
       bruker_minesweeper_strong bruker_eh bruker_background)

rm -rf "$STAGE"; mkdir -p "$STAGE"
files=()
for label in "${ORDER[@]}"; do
  src="${ARMS[$label]}"
  [ -d "$src" ] || { echo "missing input for '$label': $src" >&2; exit 1; }
  ln -s "$src" "$STAGE/$label.d"
  files+=("$STAGE/$label.d")
done
echo "staged ${#files[@]} arms under $STAGE"

mkdir -p "$OUT"
echo "running Sage LFQ over all arms (batch ${SAGE_BATCH:-3}) ..."
"$SAGE_BIN" "$CFG" -f "$FASTA" -o "$OUT" --batch-size "${SAGE_BATCH:-3}" "${files[@]}"
echo "done -> $OUT/{results.sage.tsv, lfq.tsv}"
