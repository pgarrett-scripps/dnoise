#!/usr/bin/env bash
# Re-run the diaPASEF benchmark on the FULL predicted spectral library.
#
# WHY: the default DIA path (DIANN_REUSE_RAW_LIB=1) searches the denoised/msms
# arms against the *original* arm's empirical 2nd-pass library, which scopes every
# arm to the raw-detectable precursor set -- so the denoised arms can only LOSE
# IDs, never gain. Reviewers flagged the resulting "1-4% ID cost" as an upper
# bound by construction (circular). This script instead searches each arm
# INDEPENDENTLY against the denoise-independent full predicted library
# (data/fasta/diann_lib/hybrid_diann/hybrid.predicted.speclib, ~4.84M precursors),
# so the denoised arms can gain or lose on their own merits: a fair, symmetric
# measurement rather than an upper bound.
#
# WHAT IT DOES (idempotent + reversible):
#   1. Sanity-checks the DIA-NN binary, FASTA, cached predicted lib, input .d
#      folders, and that the ORIGINAL arm's existing results were themselves
#      searched against the full predicted lib (so we keep them; raw is unchanged).
#   2. Backs up the current raw-lib-scoped denoised+msms results to
#      results/_rawlib_scoped_backup/<ds>/ (preserved as a secondary analysis),
#      then removes them so the search re-creates them.
#   3. Runs DIANN_FULL_LIB=1 11_diann.sh per DIA dataset. The original arm is
#      auto-skipped by 11_diann.sh's exists-guard (we never set DIANN_FORCE and
#      never touch original/), so only denoised + msms are re-searched.
#   4. Regenerates DIA analysis, fig4_dia, the DIA SI tables, and the DIA SI
#      figures (violins + CV).
#
# COST: full-lib searches with MBR are slow. The original arm took ~247 min
# (5-min gradient) and ~475 min (15-min) per the existing logs, so re-searching
# denoised+msms is ~8 h (5-min) + ~16 h (15-min) ~= ~24 h sequential. Run it in
# the background / overnight.
#
# AFTER IT FINISHES: the DIA numbers will have changed -- update the DIA prose in
# paper/paper.typ (abstract diaPASEF clause, Section 3.3, fig:dia caption) and
# paper/supplementary.typ (Section S5 + S1 library methods), then recompile.
#
# Usage:
#   bash scripts/52_dia_fulllib_rerun.sh                 # both gradients
#   DIA_DATASETS="dia_5min" bash scripts/52_dia_fulllib_rerun.sh   # one gradient
#   SKIP_REGEN=1 bash scripts/52_dia_fulllib_rerun.sh    # searches only, no figures
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH="$(dirname "$HERE")"
cd "$BENCH"

DIA_DATASETS="${DIA_DATASETS:-dia_5min dia_15min}"
DIANN_BIN="${DIANN_BIN:-/home/patrick-garrett/tools/diann-2.2.0/diann-linux}"
FASTA_DIANN="$BENCH/data/fasta/hybrid_diann.fasta"
LIB="$BENCH/data/fasta/diann_lib/hybrid_diann/hybrid.predicted.speclib"
BACKUP_ROOT="$BENCH/results/_rawlib_scoped_backup"

say() { echo -e "\n=== $* ==="; }

# ---------------------------------------------------------------------------
# 1. Pre-flight sanity checks
# ---------------------------------------------------------------------------
say "1. Pre-flight checks"
[ -x "$DIANN_BIN" ] || { echo "DIA-NN binary missing/not executable: $DIANN_BIN" >&2; exit 1; }
[ -f "$FASTA_DIANN" ] || { echo "stripped DIA-NN FASTA missing: $FASTA_DIANN" >&2; exit 1; }
[ -f "$LIB" ] || { echo "cached predicted lib missing: $LIB (would trigger ~8h prediction)" >&2; exit 1; }
# The predicted lib must be the FULL proteome (not a DDA-restricted allowlist).
echo "predicted lib: $LIB ($(du -h "$LIB" | cut -f1))"
head -1 "$FASTA_DIANN" | grep -q '^>sp|' \
  || { echo "FASTA does not start with >sp| headers (DIA-NN would collapse proteins)" >&2; exit 1; }

for ds in $DIA_DATASETS; do
  DATASET="$ds"; source scripts/_dataset.sh
  for v in RAW DEN MSMS; do
    dir="${!v}"
    n=$(ls -d "$dir"/*.d 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] || { echo "$ds: no .d folders in $v ($dir)" >&2; exit 1; }
    echo "$ds $v: $n .d"
  done
  # Original baseline must already be a FULL-predicted-lib search (we keep it).
  olog="$RES_ORIGINAL/diann.log"
  [ -f "$RES_ORIGINAL/report.parquet" ] || { echo "$ds: original/report.parquet missing -- cannot keep it as baseline" >&2; exit 1; }
  if ! grep -q "hybrid_diann/hybrid.predicted.speclib" "$olog" 2>/dev/null; then
    echo "$ds: original arm was NOT searched against the full predicted lib (per $olog)." >&2
    echo "      Re-run the original arm too: remove $RES_ORIGINAL/report.parquet and let this script search it." >&2
    exit 1
  fi
  echo "$ds: original baseline verified (full predicted lib) -- will keep as-is"
done

# ---------------------------------------------------------------------------
# 2. Back up + clear the raw-lib-scoped denoised/msms results
# ---------------------------------------------------------------------------
say "2. Backing up raw-lib-scoped denoised/msms results"
for ds in $DIA_DATASETS; do
  DATASET="$ds"; source scripts/_dataset.sh
  for arm in denoised msms; do
    src="$BENCH/results/$ds/$arm"
    dst="$BACKUP_ROOT/$ds/$arm"
    if [ -e "$dst" ]; then
      # Backup already taken on a prior run. NEVER touch the live dir now: on a
      # resume it holds either completed full-lib results (keep; run_arm skips
      # via its exists-guard) or a partial/absent search (run_arm re-creates it).
      echo "backup already exists ($dst); leaving live $src untouched (resume-safe)"
    elif [ -d "$src" ] && [ -f "$src/report.parquet" ]; then
      mkdir -p "$(dirname "$dst")"
      mv "$src" "$dst"
      echo "moved raw-lib results $src -> $dst"
    else
      echo "no raw-lib results to back up at $src (will be created)"
    fi
  done
done

# ---------------------------------------------------------------------------
# 3. Full-lib search (original auto-skipped via exists-guard)
# ---------------------------------------------------------------------------
say "3. DIA-NN full predicted-library search (denoised + msms)"
for ds in $DIA_DATASETS; do
  say "   $ds"
  DATASET="$ds" DIANN_FULL_LIB=1 DIANN_BIN="$DIANN_BIN" bash scripts/11_diann.sh
done

# ---------------------------------------------------------------------------
# 4. Regenerate DIA analysis, figures, and tables
# ---------------------------------------------------------------------------
if [ "${SKIP_REGEN:-0}" = "1" ]; then
  say "4. SKIP_REGEN=1 -- searches done; skipping analysis/figures"
  echo "Run later: uv run scripts/12_analyze_dia.py (per dataset), 22_paper_figures.py, 17_si_tables.py"
  exit 0
fi

say "4. Regenerating DIA analysis + figures + tables"
for ds in $DIA_DATASETS; do
  echo "--- analyze $ds ---"
  DATASET="$ds" uv run scripts/12_analyze_dia.py
done
uv run scripts/22_paper_figures.py
uv run scripts/17_si_tables.py dda_5min dda_15min dia_5min dia_15min
for g in 5min 15min; do
  ds="dia_$g"
  case " $DIA_DATASETS " in *" $ds "*)
    cp "results/$ds/analysis/lfq_ratio_violins.png" "../paper/figures/si_violins_dia_$g.png"
    cp "results/$ds/analysis/lfq_cv.png"            "../paper/figures/si_cv_dia_$g.png"
    echo "copied SI DIA figures for $g"
  ;; esac
done

say "DONE"
echo "Next: update DIA prose (paper abstract diaPASEF clause, Section 3.3, fig:dia"
echo "caption; supplementary S5 + S1) from the new numbers, then recompile both PDFs."
echo "Raw-lib-scoped results are preserved under: $BACKUP_ROOT"
