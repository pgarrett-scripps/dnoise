#!/usr/bin/env bash
# Denoise every raw .d in the dataset. Arms (select via DENOISE_ARMS):
#   ms1   -> data/<DATASET>/denoised        (vertical+halo on MS1 only)
#   msms  -> data/<DATASET>/denoised_msms   (MS1 + per-precursor MS/MS denoise)
#   wshed -> data/<DATASET>/watershed       (MS1 vertical+halo + watershed centroider, Fig 6)
# ms1/msms use config/dnoise.toml (msms adds --denoise-msms); wshed uses
# config/dnoise.watershed.toml. Default arms are "ms1 msms"; opt into the
# watershed arm with DENOISE_ARMS=wshed (or "ms1 msms wshed").
#
# STALE-FILE GUARD: each output .d carries a .dnoise_stamp = sha256(config) +
# dnoise version + arm. An output is reused only if its stamp matches the
# current config/binary/arm; otherwise it is rebuilt from scratch. This is what
# prevents the "skip (exists)" trap that silently preserved files from an older
# config or algorithm. Set DENOISE_FORCE=1 to rebuild everything regardless.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/_dataset.sh"

DNOISE="$ROOT/../target/release/dnoise"
CONFIG="${DNOISE_CONFIG:-$ROOT/config/dnoise.toml}"
WCONFIG="${DNOISE_WATERSHED_CONFIG:-$ROOT/config/dnoise.watershed.toml}"
ICONFIG="${DNOISE_INTENSITY_CONFIG:-$ROOT/config/dnoise.intensity.toml}"
DENOISE_ARMS="${DENOISE_ARMS:-ms1 msms}"   # subset of {ms1 msms wshed intensity}
FORCE="${DENOISE_FORCE:-0}"

echo "dataset: $DATASET"
echo "config : $CONFIG"
echo "arms   : $DENOISE_ARMS"
[ -f "$CONFIG" ] || { echo "config not found: $CONFIG" >&2; exit 1; }
[ -d "$RAW" ] || { echo "no raw dir: $RAW -- run 03/04_unzip first" >&2; exit 1; }

if [ ! -x "$DNOISE" ]; then
  echo "building dnoise release binary ..."
  (cd "$ROOT/.." && cargo build --release)
fi

VER="$("$DNOISE" --version)"

# denoise_arm <out_dir> <arm_label> <config> [extra dnoise flags...]
denoise_arm() {
  local outdir="$1" arm="$2" cfg="$3"; shift 3
  [ -f "$cfg" ] || { echo "config not found: $cfg" >&2; exit 1; }
  local cfg_hash; cfg_hash="$(sha256sum "$cfg" | cut -d' ' -f1)"
  local want="$VER | $cfg_hash | $arm | $*"
  mkdir -p "$outdir"
  echo "=== arm '$arm' -> $outdir ==="
  shopt -s nullglob
  for d in "$RAW"/*.d; do
    local name out stamp
    name="$(basename "$d")"
    out="$outdir/$name"
    stamp="$out/.dnoise_stamp"
    if [ "$FORCE" != "1" ] && [ -d "$out" ] && [ -f "$stamp" ] \
       && [ "$(cat "$stamp")" = "$want" ]; then
      echo "skip (up to date): $name"
      continue
    fi
    if [ -d "$out" ]; then
      echo "stale/missing stamp -- rebuilding: $name"
      rm -rf "$out"
    fi
    echo "denoising $name ..."
    "$DNOISE" "$d" "$out" --config "$cfg" "$@"
    printf '%s' "$want" > "$stamp"
  done
}

# Acquisition-aware noise gates, applied at denoise time (benchmark default):
#   DDA -> drop MS1 outside the ddaPASEF/PASEF selection polygon (--ms1-polygon).
#   DIA -> drop MS1 (--dia-ms1-window) and, on the MS/MS arm, MS/MS
#          (--dia-window) points outside the diaPASEF isolation windows.
# The polygon is DDA-only: diaPASEF stores a very restrictive polygon (~91% of
# survey MS1) that would starve DIA-NN's MS1 quant, so DIA uses the per-window
# gates instead. All gates are safe no-ops when their tables/polygon are absent.
# --dia-window goes on the msms arm only, so the MS1-only "denoised" arm stays
# MS1-only. Gate flags enter the .dnoise_stamp, so arms rebuild when they change.
case "$DATASET" in
  dia*) MS1_GATE=(--dia-ms1-window); MSMS_GATE=(--dia-ms1-window --dia-window) ;;
  *)    MS1_GATE=(--ms1-polygon);    MSMS_GATE=(--ms1-polygon) ;;
esac

for arm in $DENOISE_ARMS; do
  case "$arm" in
    ms1)       denoise_arm "$DEN"   ms1       "$CONFIG"  "${MS1_GATE[@]}" ;;
    msms)      denoise_arm "$MSMS"  msms      "$CONFIG"  --denoise-msms "${MSMS_GATE[@]}" ;;
    wshed)     denoise_arm "$WSHED" wshed     "$WCONFIG" "${MS1_GATE[@]}" ;;
    intensity) denoise_arm "$INT"   intensity "$ICONFIG" "${MS1_GATE[@]}" ;;
    *) echo "unknown arm: $arm" >&2; exit 1 ;;
  esac
done

echo "done."
