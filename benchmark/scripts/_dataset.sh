# Shared dataset namespace for the benchmark scripts. Source this, don't run it.
#
#   DATASET selects which acquisition/gradient set to operate on. Everything is
#   namespaced under data/<DATASET>/ and results/<DATASET>/ so multiple datasets
#   (dda_5min, dda_15min, dia_5min, ...) coexist without clobbering each other.
#   Default is dda_5min for backwards compatibility.
#
#   config/files.<DATASET>.tsv  -- the file list for the dataset
#   config/dnoise.toml          -- SHARED denoise params (gradient-independent)
#   config/sage.json            -- SHARED Sage config (FASTA + tolerances)
#   data/fasta/hybrid.fasta     -- SHARED hybrid proteome
#
# Callers get: $ROOT $DATASET $FILES $ZIPS $RAW $DEN $MSMS $RES_* $FASTA
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATASET="${DATASET:-dda_5min}"

FILES="$ROOT/config/files.$DATASET.tsv"
# FASTA is overridable per dataset (e.g. the UPS2+E.coli DB for the Pro 2 arm);
# defaults to the shared HYE hybrid proteome.
FASTA="${FASTA:-$ROOT/data/fasta/hybrid.fasta}"

DATA="$ROOT/data/$DATASET"
ZIPS="$DATA/zips"
RAW="$DATA/raw"
DEN="$DATA/denoised"          # MS1-only denoise arm (streak: vertical+halo)
MSMS="$DATA/denoised_msms"    # MS1 + MS/MS denoise arm
WSHED="$DATA/watershed"       # MS1 vertical+halo + watershed centroider (Fig 6)
INT="$DATA/denoised_intensity" # MS1 strict-intensity-threshold control arm
INTMSMS="$DATA/denoised_intensity_msms" # MS1+MS/MS strict-intensity-threshold control arm

RESULTS="$ROOT/results/$DATASET"
RES_ORIGINAL="$RESULTS/original"
RES_DENOISED="$RESULTS/denoised"
RES_MSMS="$RESULTS/msms"
RES_WSHED="$RESULTS/watershed"
RES_INTENSITY="$RESULTS/intensity"
RES_INTMSMS="$RESULTS/intensity_msms"
RES_ANALYSIS="$RESULTS/analysis"

export DATASET   # so child python (06/07) picks up the same namespace
