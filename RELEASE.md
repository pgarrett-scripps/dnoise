# Releasing v0.5.0

Status: **in progress, 2026-09-25.** v0.4.0 was released 2026-09-23 (9d2054c).

From 0.5.0 the workspace publishes two crates, `dnoise-core` and `dnoise`, at
the same version. Bump `dnoise-core/Cargo.toml`'s `version` and the
`dnoise-core` dependency's `version` in `Cargo.toml` with the rest, and run
`cargo publish --dry-run --locked -p dnoise-core` in the local gate.
`release.yml` publishes `dnoise-core` before `dnoise`.

## Steps

1. Version 0.5.0 in `Cargo.toml`, `dnoise-core/Cargo.toml`, `dnoise-gui/Cargo.toml` and `Cargo.lock`,
   then `just cite-sync` to copy it and today's date into `CITATION.cff`;
   `CHANGELOG.md` has a dated 0.5.0 section and `[Unreleased]` is empty.
2. Local gate: `cargo fmt --all --check`, `cargo clippy --workspace
   --all-targets -- -D warnings`, `cargo test --release --workspace`,
   `cargo publish --dry-run --locked`.
3. Push `release/0.5.0`, open a PR to `main`, wait for CI to pass, then
   fast-forward `main` so the tagged commit is the one tested.
4. `gh release create v0.5.0 --target main` with the 0.5.0 changelog as notes.
   `release.yml` attaches the binaries and publishes `dnoise-core` then `dnoise`; Zenodo mints the
   version DOI.
5. Cite the version DOI in the paper. `CITATION.cff` and `README.md` use the
   concept DOI, so no DOI changes there.

## Output changes from 0.3.0

Default MS1 output changes on both acquisition types (feature-level overlap
gates, no fixed pads, Bruker-calibrated 1/K0). Dry-run on the 5-minute benchmark
runs, all points kept: ddaPASEF 111,722,871 of 297,755,895; diaPASEF
755,841,088 of 1,321,535,492. `--linear-mobility` reproduces the
pre-calibration counts. MS/MS output is unchanged.

---

# v0.3.0 release record

Status: **released on 2026-09-14.** GitHub binaries, the crates.io package,
and the Zenodo archive are published.

- [GitHub v0.3.0](https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.3.0)
  tags commit `f14fbc22f301f9600b93da62560bac6bbea19b04`.
- [Zenodo v0.3.0](https://doi.org/10.5281/zenodo.22756840) is published; the
  downloaded archive's checksum, source commit, version, and Windows fix were verified.
- The all-version DOI remains `10.5281/zenodo.21959649`. The immutable v0.1.0
  archive is `10.5281/zenodo.21959650`; use that record for the original paper results.
- [crates.io v0.3.0](https://crates.io/crates/dnoise/0.3.0) was published by
  [recovery run 34887161343](https://github.com/pgarrett-scripps/dnoise/actions/runs/34887161343).

## Still to do

1. **Re-derive the paper's diaPASEF MS/MS arm.** 0.3.0 fixes a v0.1.0 defect in
   that arm (see below), so the manuscript's optional-fragment-denoising numbers
   were produced by the buggy path. The paper currently says so in its
   availability statement and SI. Tracked in the paper repo's `TODO.md`, which
   has the ordered steps. This does not block the software release, only the
   claim that the paper and the release agree everywhere.

## Done

- Fixed Windows batch preflight, which queried an incomplete canonical drive
  prefix. Added regression coverage for canonical paths and missing children.
  All eight [CI checks](https://github.com/pgarrett-scripps/dnoise/actions/runs/34884608764)
  passed, including formatting, Clippy, and workspace tests on Linux, macOS,
  and Windows; MSRV, library-only, stable, dependency audit, and SDK-comparison checks.
- Repeated local formatting, Clippy, and workspace tests successfully. The final
  `cargo publish --dry-run --locked` verified 89 files (949.3 KiB uncompressed).
- Published the GitHub release with the 0.3.0 changelog as its notes. Updated
  `README.md` and `CITATION.cff` with the version-specific Zenodo DOI.
- All three platform test/build jobs passed and the CLI/GUI archives are attached.
  The downloaded Linux CLI reports `dnoise 0.3.0`; its archive SHA-256 matches
  the digest on the published GitHub asset.
- Installed the published crate with `cargo install dnoise --version 0.3.0 --locked`
  into an isolated temporary directory; the installed command reports `dnoise 0.3.0`.
- Fixed a packaging-workflow issue: downloaded binary archives in `dist/` made
  Cargo reject the source checkout as dirty. Future releases stage downloads
  outside the checkout. Published this crate from the unchanged v0.3.0 tag using
  the recovery workflow, which checks the release downloads and package version.
- Updated the local paper and SI release links, preserving the v0.1.0 benchmark
  reference. The manuscript's existing unresolved diaPASEF rerun note blocks
  `just paper`; `just verify` also reports existing declaration and stale-output
  issues. These do not block the software release.
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
