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
