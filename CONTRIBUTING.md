# Contributing to dnoise

Thanks for your interest in `dnoise`. It is a Rust CLI and library that denoises
Bruker timsTOF `.d` folders with an iterative vertical ion-mobility feature
filter. Contributions of all kinds are welcome, from bug reports and
documentation fixes to new features.

## Reporting bugs and requesting features

Please use [GitHub issues](https://github.com/pgarrett-scripps/dnoise/issues).

- For bugs, include what you ran, what you expected, what happened, and enough
  detail to reproduce (OS, Rust version, a minimal `.d` sample or its
  characteristics if you can share them).
- For feature requests, describe the use case and the outcome you want.

Search existing issues first to avoid duplicates.

## Development setup

```sh
git clone https://github.com/pgarrett-scripps/dnoise.git
cd dnoise
cargo build
cargo test --workspace
```

`dnoise` targets a minimum supported Rust version (MSRV) of **1.85** and uses
**edition 2024**. Install a matching toolchain with `rustup` if needed.

## Code style

Before submitting, run:

```sh
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
```

CI enforces both. Formatting must match `cargo fmt`, and Clippy runs with
warnings treated as errors, so any warning will fail the build. Please fix or
explicitly address every lint rather than allowing it away without reason.

## Adding tests

- Put unit tests inline in the module they cover, inside a `#[cfg(test)]`
  module.
- Put integration tests in the top-level `tests/` directory.

New behavior should come with tests. Bug fixes should ideally include a test
that fails before the fix and passes after.

## Pull request process

1. Branch from `main`.
2. Keep commits focused and the history readable. One logical change per commit
   where practical.
3. Make sure `cargo build`, `cargo test --workspace`, `cargo fmt --all --check`,
   and `cargo clippy --workspace --all-targets -- -D warnings` all pass locally.
4. Open the PR against `main` and make sure CI is green.
5. Reference any related issue in the PR description.

## Algorithm reference

For the specification of the denoising algorithm itself, see
[ALGORITHM.md](ALGORITHM.md).

## License and Code of Conduct

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE), and that you will follow the
[Code of Conduct](CODE_OF_CONDUCT.md).
