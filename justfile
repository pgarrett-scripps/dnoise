# Default target
default:
  @just --list

# Build (debug)
build:
  cargo build

# Build optimized release binary
release:
  cargo build --release

# Run tests
test:
  cargo test

# Run clippy on all targets, warnings as errors
lint:
  cargo clippy --all-targets -- -D warnings

# Format code
format:
  cargo fmt

# Fail if dnoise-core gained an I/O or native dependency
core-deps:
  .github/scripts/check-core-deps.sh

# Check formatting without modifying files
format-check:
  cargo fmt --check

# Lint + format check + tests
check:
  just lint
  just format-check
  just test

# Denoise a .d folder: just denoise INPUT.d OUTPUT.d [args...]
denoise input output *args:
  cargo run --release -- "{{input}}" "{{output}}" {{args}}

# Validate a (denoised) .d folder re-reads and matches its DB
validate path:
  cargo run --release --example validate -- "{{path}}"

# Remove build artifacts
clean:
  cargo clean

# Copy Cargo.toml's version and a release date (default today) into CITATION.cff
cite-sync date=`date +%F`:
  #!/usr/bin/env bash
  set -euo pipefail
  v=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)
  sed -i -E "s/^version: .*/version: $v/; s/^date-released: .*/date-released: \"{{date}}\"/" CITATION.cff
  echo "CITATION.cff: version $v, date-released {{date}}"

# Fail if CITATION.cff's version differs from Cargo.toml's (CI runs the same check)
cite-check:
  #!/usr/bin/env bash
  set -euo pipefail
  crate=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)
  cff=$(sed -n 's/^version: //p' CITATION.cff)
  grep -Eq '^date-released: "[0-9]{4}-[0-9]{2}-[0-9]{2}"$' CITATION.cff || { echo "CITATION.cff: missing or malformed date-released"; exit 1; }
  [ "$crate" = "$cff" ] || { echo "CITATION.cff version $cff != Cargo.toml $crate: run just cite-sync"; exit 1; }
  echo "CITATION.cff matches Cargo.toml ($crate)"
