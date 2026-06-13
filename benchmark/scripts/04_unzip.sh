#!/usr/bin/env bash
# Unzip the timsTOF zips downloaded by scripts/03_download.sh into data/raw.
# Resumable (skips an existing .d). Keeps the zip files.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/_dataset.sh"
echo "dataset: $DATASET  ($FILES)"
mkdir -p "$RAW"

[ -f "$FILES" ] || { echo "missing $FILES — run scripts/01_list_files.py first"; exit 1; }

tail -n +2 "$FILES" | while IFS=$'\t' read -r cond rep name size url; do
  dpath="$RAW/${name%.zip}"
  if [ -d "$dpath" ]; then echo "skip unzip (exists): ${name%.zip}"; continue; fi
  zip="$ZIPS/$name"
  [ -f "$zip" ] || { echo "missing zip: $name — run scripts/03_download.sh first"; exit 1; }
  echo "unzipping $name ..."
  unzip -q -o "$zip" -d "$RAW"
done

echo "raw .d folders:"
ls -d "$RAW"/*.d
