#!/usr/bin/env bash
# Download the timsTOF zips listed in config/files.$DATASET.tsv (DATASET env).
#
# Connection count is configurable via CONN (default 1 = single-stream).
#
#   History: this defaulted to single-stream (-x1 -s1) because multi-connection
#   segmented transfers (aria2 -x16) were suspected of stitching byte ranges at
#   slightly wrong offsets -- file ends up the RIGHT SIZE but with corrupt bytes
#   mid-file, which a size check misses. We later tested -x16 -s16 against the
#   PRIDE HTTPS endpoint (6/6 downloads of a 5min DDA file passed `unzip -t`),
#   and saw a ~10x speedup with no corruption. So multi-connection is now opt-in
#   via CONN, e.g.  CONN=16 DATASET=dia_5min ./scripts/03_download.sh
#
#   The CRC safety net below makes this safe regardless: every zip is verified
#   with `unzip -t` (CRC of every member) before it counts as done. Corrupt or
#   incomplete zips are deleted and re-downloaded from scratch -- so even a rare
#   mis-stitched file cannot pass silently; it just gets refetched.
#
# Resumable: aria2 resumes from a dropped connection (-c). A zip that exists but
# FAILS verification is NOT skipped -- it is re-downloaded. Only a zip that
# passes `unzip -t` is treated as complete.
#
# Does NOT unzip the data -- see scripts/04_unzip.sh. Zips are kept.
set -uo pipefail   # NOTE: no -e: one failed file must not abort the whole run

source "$(dirname "${BASH_SOURCE[0]}")/_dataset.sh"
CONN="${CONN:-1}"      # aria2 connections/splits per file; CONN=16 for ~10x speed
echo "dataset: $DATASET  ($FILES)  [connections: $CONN]"
mkdir -p "$ZIPS" "$RAW"

MAX_ATTEMPTS=4          # download+verify attempts per file before giving up

[ -f "$FILES" ] || { echo "missing $FILES -- run scripts/01_list_files.py first"; exit 1; }

# verify_zip <path> <expected_size> : 0 = good, 1 = bad/missing
verify_zip() {
  local zip="$1" expected="$2" actual
  [ -f "$zip" ] || return 1
  actual="$(stat -c%s "$zip" 2>/dev/null || echo 0)"
  if [ -n "$expected" ] && [ "$expected" != "0" ] && [ "$actual" != "$expected" ]; then
    echo "  size mismatch: got $actual, expected $expected"
    return 1
  fi
  # CRC-check every member. This is what catches a right-size/wrong-bytes file.
  if ! unzip -t -q "$zip" >/dev/null 2>&1; then
    echo "  zip integrity check (unzip -t) FAILED"
    return 1
  fi
  return 0
}

failed=()

tail -n +2 "$FILES" | while IFS=$'\t' read -r cond rep name size url; do
  zip="$ZIPS/$name"
  https="${url/ftp:\/\//https://}"

  # Already have a verified-good zip? Skip it.
  if [ -f "$zip" ]; then
    echo "checking existing $name ..."
    if verify_zip "$zip" "$size"; then
      echo "skip download (verified good): $name"
      continue
    fi
    echo "existing $name is bad -- removing and re-downloading"
    rm -f "$zip" "$zip.aria2"
  fi

  ok=0
  for attempt in $(seq 1 "$MAX_ATTEMPTS"); do
    echo "downloading $name (attempt $attempt/$MAX_ATTEMPTS) ..."
    # -x/-s CONN: connections + splits per file (CONN=1 -> single-stream).
    # -k1M: min split size, so CONN connections actually parallelize big files.
    # -c: resume a dropped connection.  --max-tries/--retry-wait: ride out blips.
    aria2c -x"$CONN" -s"$CONN" -k1M -c \
      --max-tries=10 --retry-wait=30 --timeout=120 \
      --console-log-level=warn --summary-interval=30 \
      --allow-overwrite=true -d "$ZIPS" -o "$name" "$https"

    echo "verifying $name ..."
    if verify_zip "$zip" "$size"; then
      echo "OK: $name"
      ok=1
      break
    fi

    echo "verification failed for $name -- deleting and retrying from scratch"
    rm -f "$zip" "$zip.aria2"
  done

  if [ "$ok" != "1" ]; then
    echo "GIVING UP on $name after $MAX_ATTEMPTS attempts"
    failed+=("$name")
  fi
done

echo
echo "zips:"
ls -1 "$ZIPS"/*.zip 2>/dev/null || echo "(none)"

# NOTE: 'failed' is populated inside a 'while' in a pipeline (subshell), so it is
# not visible here. Do a final authoritative pass over all files instead.
echo
echo "final verification:"
bad=0
tail -n +2 "$FILES" | while IFS=$'\t' read -r cond rep name size url; do
  if verify_zip "$ZIPS/$name" "$size" >/dev/null 2>&1; then
    echo "  OK   $name"
  else
    echo "  BAD  $name"
  fi
done
