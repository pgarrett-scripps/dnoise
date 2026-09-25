#!/usr/bin/env bash
# Fail if dnoise-core's normal (non-dev, non-build) dependency tree pulls in
# I/O or native code: timsrust*, rusqlite, libsqlite3-sys, zstd*, or any crate
# with a `links` key (rayon-core has one, so rayon is out too). sage-plus and
# other embedders rely on this.
set -euo pipefail
cd "$(dirname "$0")/../.."

tree=$(cargo tree -p dnoise-core -e normal --prefix none --locked)
echo "$tree"

forbidden='^(timsrust[a-z0-9_-]*|rusqlite|libsqlite3-sys|zstd|zstd-sys|zstd-safe|rayon|rayon-core) '
if echo "$tree" | grep -Eq "$forbidden"; then
  echo "error: forbidden dependency in dnoise-core:" >&2
  echo "$tree" | grep -E "$forbidden" >&2
  exit 1
fi

# Packages in the tree (name + version) that declare a `links` key.
echo "$tree" | awk '{print $1, $2}' | sort -u > "${TMPDIR:-/tmp}/dnoise-core-deps.txt"
cargo metadata --format-version 1 --locked | python3 -c '
import json, sys
deps = {tuple(l.split()) for l in open(sys.argv[1]) if l.strip()}
meta = json.load(sys.stdin)
bad = sorted(p["name"] for p in meta["packages"]
             if (p["name"], "v" + p["version"]) in deps and p.get("links"))
if bad:
    sys.exit("error: dnoise-core depends on crates with a links key: " + ", ".join(bad))
print(f"dnoise-core: {len(deps)} normal dependencies, none forbidden, none with links")
' "${TMPDIR:-/tmp}/dnoise-core-deps.txt"
