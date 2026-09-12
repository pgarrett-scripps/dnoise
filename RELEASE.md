# Releasing v0.3.0

Status: **not released.** Everything below the line is done; everything above it
is not. Nothing here has been tagged, published, or archived.

## Still to do

1. **Re-derive the paper's diaPASEF MS/MS arm.** 0.3.0 fixes a v0.1.0 defect in
   that arm (see below), so the manuscript's optional-fragment-denoising numbers
   were produced by the buggy path. The paper currently says so in its
   availability statement and SI. Tracked in the paper repo's `TODO.md`, which
   has the ordered steps. This does not block the software release, only the
   claim that the paper and the release agree everywhere.
2. Publish a GitHub release for `v0.3.0`, with the 0.3.0 section of
   CHANGELOG.md as its notes. Creating the release creates the tag, so there is
   nothing to tag first. That fires `.github/workflows/release.yml`, which
   builds the per-OS binaries, attaches them, and publishes the crate. Needs the
   `CARGO_REGISTRY_TOKEN` repository secret. A draft release does not fire it;
   publishing the draft later does.
3. Confirm the Zenodo webhook minted a record for the tag, then put its DOI in:
   - `README.md` (citation block, marked with a `RELEASE:` comment)
   - `CITATION.cff` (`doi:` under `identifiers`)
   - the paper repo's `paper.typ` and `si-body.typ`, at the two
     `FIXME(release)` comments
4. Verify the published crate installs clean: `cargo install dnoise --version 0.3.0`.

## Done

- Version bumped to 0.3.0 in `Cargo.toml`, `dnoise-gui/Cargo.toml`, `CITATION.cff`
  (with `date-released`).
- `CHANGELOG.md` has a dated 0.3.0 section; `[Unreleased]` is empty.
- `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`,
  and `cargo test --release` (161 tests) all pass.
- `cargo publish --dry-run` packages cleanly (87 files, 939 KiB).
- Equivalence against v0.1.0 checked on real acquisitions, as CONTRIBUTING
  requires before a writer change. Built v0.1.0 from its tag and ran both
  binaries over the benchmark with the paper's configs, comparing
  `analysis.tdf` and `analysis.tdf_bin`:

  | Arm | Pairs | Result |
  |---|---|---|
  | ddaPASEF MS1-only / MS1+MS/MS / watershed, both gradients | 24 | byte-identical |
  | diaPASEF MS1-only | 6 | byte-identical |
  | diaPASEF MS1+MS/MS | 6 | differs, deliberately |

  The diaPASEF MS/MS difference is the v0.1.0 fix: `filter_per_window` is meant
  to filter inside each isolation window separately, since neighbouring windows
  isolate unrelated precursor m/z bands, but v0.1.0's `read_dia_windows`
  coalesced windows that touch in scan space. The benchmark's windows touch
  exactly, so all three per group merged into one interval and the filter ran
  over the whole frame. 0.3.0 keeps static windows separate and retains about
  0.37% fewer MS/MS points. MS1 is unaffected.
- Empty-frame encoding fixed and covered by a regression test. 0.2.0-dev wrote
  0-peak frames as a header-only 8-byte record that timsrust cannot decode,
  which made any file with an emptied frame unreadable to Sage and every other
  timsrust consumer.
