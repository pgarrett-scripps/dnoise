# Releasing dnoise

**No release is in progress.** The latest release is v0.6.0, tagged 2026-10-11 at
`3c1772e`, archived as
[10.5281/zenodo.23301791](https://doi.org/10.5281/zenodo.23301791).

The workspace publishes two crates at the same version, `dnoise-core` and
`dnoise`, and `release.yml` publishes `dnoise-core` first. `dnoise-gui` is not
published to crates.io and ships only as an attached binary archive.

## Steps for version X.Y.Z

1. Set X.Y.Z in `Cargo.toml`, `dnoise-core/Cargo.toml`, `dnoise-gui/Cargo.toml`,
   the `dnoise-core` dependency's `version` in `Cargo.toml`, and `Cargo.lock`.
   Run `just cite-sync` to copy the version and today's date into `CITATION.cff`,
   which CI's citation job checks against `Cargo.toml`. `CHANGELOG.md` needs a
   dated X.Y.Z section and an empty `[Unreleased]`.
2. Local gate: `cargo fmt --all --check`, `cargo clippy --workspace
   --all-targets -- -D warnings`, `cargo test --release --workspace`,
   `cargo publish --dry-run --locked -p dnoise-core`, then
   `cargo publish --dry-run --locked --workspace`.
3. Push `release/X.Y.Z`, open a PR to `main`, wait for CI to pass, then
   fast-forward `main` so the tagged commit is the one CI tested.
4. `gh release create vX.Y.Z --target main` with that version's changelog section
   as the notes. `release.yml` attaches the binaries, publishes `dnoise-core` and
   then `dnoise`, and Zenodo mints the version DOI. A pending DOI or a 409 is
   normal, so never release the same version twice. Stage any downloaded binary
   archives outside the checkout: a populated `dist/` once made Cargo reject the
   source as dirty.
5. Add the version DOI to the table below and cite it from whatever is pinned to
   that version. `CITATION.cff` and `README.md` carry the concept DOI
   `10.5281/zenodo.21959649`, which resolves to the latest release, so they need
   no DOI change per release.

## Released versions

| Version | Date | Tag commit | Version DOI | Default output |
|---|---|---|---|---|
| 0.6.0 | 2026-10-11 | `3c1772e` | [10.5281/zenodo.23301791](https://doi.org/10.5281/zenodo.23301791) | unchanged from 0.5.0 |
| 0.5.0 | 2026-09-25 | `9beeb3d` | [10.5281/zenodo.22965157](https://doi.org/10.5281/zenodo.22965157) | **changed** on both acquisition types |
| 0.4.0 | 2026-09-23 | `9d2054c` | [10.5281/zenodo.22929326](https://doi.org/10.5281/zenodo.22929326) | **changed** on both acquisition types |
| 0.3.0 | 2026-09-14 | `f14fbc2` | [10.5281/zenodo.22756840](https://doi.org/10.5281/zenodo.22756840) | **changed** on diaPASEF MS/MS |
| 0.1.0 | 2026-08-15 | `fa750c9` | [10.5281/zenodo.21959650](https://doi.org/10.5281/zenodo.21959650) | original paper submission |

The manuscript is pinned to v0.4.0. Use the v0.1.0 record for the original
submission's results. Each "changed" row is spelled out in that version's
[CHANGELOG.md](CHANGELOG.md) entry, with the flags that restore the previous
defaults where they can be restored.

Within 0.5.0 and 0.6.0, two changes were checked to leave output alone: the
`dnoise-core` crate split and the move off `timsrust-tdf` are byte-identical to
0.4.0 on the 5-minute benchmark runs. What changed the output in 0.5.0 was the
halo `peak_fraction` default of 0.10, the MS1 gate pads of 3 Th and 0.015 1/K0
with unlimited feature reach, and the MS1 stage order of streak, gate, halo.

---

# Release records

## v0.3.0

Released on 2026-09-14. GitHub binaries, the crates.io package, and the Zenodo
archive are published.

- [GitHub v0.3.0](https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.3.0)
  tags commit `f14fbc22f301f9600b93da62560bac6bbea19b04`.
- [Zenodo v0.3.0](https://doi.org/10.5281/zenodo.22756840) is published. The
  downloaded archive's checksum, source commit, version, and Windows fix were verified.
- [crates.io v0.3.0](https://crates.io/crates/dnoise/0.3.0) was published by
  [recovery run 34887161343](https://github.com/pgarrett-scripps/dnoise/actions/runs/34887161343).
- The diaPASEF MS/MS arm of the paper was produced with v0.1.0, before the
  `read_dia_windows` fix below. Decision of 2026-09-22: that arm is not rerun,
  and the manuscript discloses it in the text instead. Nothing here blocks a
  software release.

### Verified at the time

- Fixed Windows batch preflight, which queried an incomplete canonical drive
  prefix. Added regression coverage for canonical paths and missing children.
  All eight [CI checks](https://github.com/pgarrett-scripps/dnoise/actions/runs/34884608764)
  passed, including formatting, Clippy, and workspace tests on Linux, macOS,
  and Windows, plus MSRV, library-only, stable, dependency audit, and
  SDK-comparison checks.
- Repeated local formatting, Clippy, and workspace tests successfully. The final
  `cargo publish --dry-run --locked` verified 89 files (949.3 KiB uncompressed).
- Published the GitHub release with the 0.3.0 changelog as its notes.
- All three platform test/build jobs passed and the CLI/GUI archives are attached.
  The downloaded Linux CLI reports `dnoise 0.3.0`, and its archive SHA-256 matches
  the digest on the published GitHub asset.
- Installed the published crate with `cargo install dnoise --version 0.3.0 --locked`
  into an isolated temporary directory. The installed command reports `dnoise 0.3.0`.
- Fixed a packaging-workflow issue: downloaded binary archives in `dist/` made
  Cargo reject the source checkout as dirty. Later releases stage downloads
  outside the checkout. Published this crate from the unchanged v0.3.0 tag using
  the recovery workflow, which checks the release downloads and package version.
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
