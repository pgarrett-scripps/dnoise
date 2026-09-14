# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2026-09-14

Adds prm-PASEF and the scanning diaPASEF variants, a streaming library API, and
in-file provenance.

Denoised output is byte-identical to 0.1.0 for ddaPASEF on all three arms
(MS1-only, MS1 + MS/MS, watershed centroider) and for diaPASEF MS1-only,
verified by running both builds over the benchmark acquisitions and comparing
the stored frame data and metadata.

**diaPASEF MS1 + MS/MS output changed, because 0.1.0 was wrong here.**
`filter_per_window` exists to filter fragments inside each isolation window
separately, since neighbouring windows isolate unrelated precursor m/z bands.
0.1.0's `read_dia_windows` coalesced any window whose `ScanNumBegin` was at or
before the previous `ScanNumEnd`, and windows that merely touch satisfy that.
Where a group's windows touch, they all merged into one interval and the filter
ran over the whole frame, letting a fragment run span a quadrupole jump between
unrelated m/z bands. 0.3.0 keeps static windows separate and merges only a
genuine scanning ramp: consecutive one-scan steps of equal width moving
monotonically in m/z. About 0.37% fewer MS/MS points are retained, across
roughly half the MS/MS frames, and the points now removed are mostly noise runs
that used to bridge a window boundary. MS1 is unaffected. Anyone reproducing
published diaPASEF fragment-denoising numbers from 0.1.0 should pin 0.1.0 or
re-derive them.

### Fixed
- Resolve canonical Windows paths without querying an incomplete drive prefix,
  fixing batch manifests that failed with "Incorrect function".
- Write empty (0-peak) frames the long way, as a header plus a compressed
  all-zero scan table, rather than the header-only 8-byte record. timsrust
  zstd-decodes whatever follows the header unconditionally and fails on an
  absent payload, so a single emptied frame made the whole `.d` unreadable to
  every downstream consumer, Sage included. MS/MS denoising empties frames
  routinely. The reader still accepts both shapes.

### Added and changed
- Reject physical gates/crops that require unsupported multi-calibration conversions;
  retain raw-coordinate filtering and RT/intensity-only crops.
- Report separate MS1/MS/MS intensity retention and actual neighbor usage, including
  intensity removed by RT cropping in the input totals.
- Add a reusable real-acquisition comparison command with independent SDK checks,
  pinned input checksums, exact executable hashes, recipes, and failure reports.

- Add prm-PASEF support: validate PRM target/event metadata, detect PRM and mixed
  acquisitions in CLI/GUI/reports, and preserve fragment points by default.
  `--denoise-msms` enables experimental filtering within each PRM isolation
  event, including event-local halo/smoothing/centroiding. Cross-frame pooling
  is disabled by default. Record the active path in provenance. Disable discovery gates on
  PRM/mixed/unknown data; mixed/unknown runs reject fragment denoising. Validate a
  public PRM run with the Bruker SDK; see [targeted support](docs/targeted.md).
- Add experimental PRM/DIA/MS1 neighbor radii in CLI, TOML, library stages, GUI,
  and provenance. Sum compatible nearby observations for keep decisions without
  importing peaks or intensities. Enforce event boundaries, calibration segments,
  acquisition settings, and a configurable RT-distance limit (default 5 seconds).
  Keep `--frame-half-width` / `frame_half_width` as MS1 aliases; correct obsolete
  averaging documentation. See [neighbor support](docs/neighbors.md).
- Support public synchro-PASEF, midia-PASEF, and Slice-PASEF examples through the
  DIA path. Keep static touching windows separate, group continuous one-scan
  trajectories for filtering, and match complete trajectories for temporal
  support. Validate DIA references before opening the decoder. Record scanning
  filtering with `ActiveGates::dia_scan_varying`; see [examples](docs/acquisition-examples.md).
- Protect input/output paths and install only validated temporary outputs, with
  ownership-aware cancellation cleanup and recoverable in-place replacement.
- Store processing history and exact configuration inside every completed `.d`.
  Add shared configuration resolution, versioned reports, checked codec APIs,
  synthetic DDA/DIA integration tests, and `validate`/`metadata` commands.
- Add native folder pickers, batch collision checks, retry controls, completed
  frame comparisons, and per-frame cancellation checks. Expose frame batch size.
- Add portable JSON batch manifests and `dnoise batch`, including per-file results.
- Upgrade affected dependencies and the GUI framework. Add dependency audits,
  monthly update PRs, explicit compiler selection, and platform test gates.
- Parallelize full-frame validation with Rayon using bounded read batches and
  at most four decoder tasks; retain structural, content, and cancellation checks.
- Make full input/output validation optional via `--skip-validation`, shared
  configuration, and a GUI checkbox. Keep full validation on by default, retain
  structural/file-safety checks, and record the mode in provenance and recipes.

Library migration: `Stages` adds `neighbors: NeighborParams`; `ActiveGates` adds
`ms1_neighbors`, `prm_neighbors`, and `dia_neighbors`. `RunOptions` adds
`frame_batch_size` and `skip_validation`; use
`..RunOptions::default()` when constructing options. The canonical empty codec
record now decodes successfully. Reports use schema version 1. The GUI requires
Rust 1.95+, while the root CLI/library minimum remains 1.85.
`Acquisition` adds `PrmPasef` and `Mixed`; update exhaustive matches. Detection
now validates PRM metadata and classifies multiple nonzero frame types as mixed
instead of choosing DDA or DIA by precedence.
`ActiveGates` adds `prm_per_event`; use its default when constructing partial
values. Reports retain schema version 1 with this additive field.


## [0.1.0] - 2026-08-14

First public release.

### Added
- **Region-of-interest crop / trim** (`crop` module, `CropParams`, `CropGate`):
  keep only points inside an axis-aligned `(m/z, 1/K0, intensity)` box plus a
  retention-time window, to carve a smaller `.d` out of a large one. CLI:
  `--mz-min`/`--mz-max` (Da), `--im-min`/`--im-max` (1/K0),
  `--min-intensity`/`--max-intensity`, and `--rt-min`/`--rt-max` (minutes). m/z and
  mobility become integer `(TOF, scan)` ranges via the run calibration and apply to
  every frame. RT bounds empty out-of-window frames (never deleted, so the frame
  axis and dependent tables stay valid and SDK-compatible). Composes with the
  denoiser, or run `--crop-only` to apply just the crop with no denoising.
- **Acquisition-aware noise gates, on by default** (`Ms1PolygonParams`,
  `DiaMs1WindowParams`, `DdaWindowParams`, `Acquisition`, `detect_acquisition`).
  Each gate auto-detects its defining geometry from the run and is a silent
  no-op when it is absent, so a single default set picks the right gate per
  acquisition. MS1 gates, always on by default: `--ms1-polygon`
  (ddaPASEF/PASEF) drops MS1 points outside the run's precursor-selection
  polygon. `--dia-ms1-window` (diaPASEF) drops MS1 points whose
  `(m/z, mobility)` falls outside every isolation window. Both pad their
  boundary in physical units (default 5 Da / 0.05 1/K0) so an edge precursor
  keeps its isotopic envelope and mobility spread. MS/MS gates, on by default
  whenever MS/MS frames are filtered (`--denoise-msms` / `--all-frames`):
  `--dda-window` (ddaPASEF) and `--dia-window` (diaPASEF) drop MS/MS points
  whose mobility scan falls outside every isolation event/window for their
  frame (`--dda-window-scan-pad` / `--dia-window-scan-pad` add edge leniency).
  On standard ddaPASEF files `--dda-window` removes nothing and just enforces
  the invariant. `--no-*` flags (or `<key> = false` in the config) force any
  gate off. An explicit flag or config value always wins over the default.
  `Stages` carries the per-gate `Option<&…Params>` fields.
- **Calibration-segment warning** (`tdf::count_calibration_segments`): when an
  acquisition gate is active and the run's frames reference more than one
  `MzCalibration` or `TimsCalibration` segment, dnoise logs a warning (in both
  the file writer and `RunContext`). The gates convert their physical-unit
  definitions to index space once with the run-level calibration timsrust
  exposes, which is exact only for a single-segment run. On a multi-segment
  file the gate boundary would be ppm-scale offset on frames referencing other
  segments, and the warning says so instead of gating silently.
- **Pinned toolchain** (`rust-toolchain.toml`, 1.97.1): the exact toolchain the
  paper binary was built with. `rust-version = "1.85"` remains the MSRV.
- **Dry runs, sampling & JSON reports**: `--dry-run` runs the full pipeline and
  reports the reduction without writing any output. `--sample <f>` (dry-run only)
  processes a deterministic fraction of frames for a fast estimate
  (`--sample-seed`). `--report <FILE>` writes the effective config plus per-MS-level
  reduction statistics as JSON. Library: `denoise_with_options`, `RunOptions`,
  `SampleSpec`, and a widened `DenoiseStats` (per-level point counts, summed
  intensities, cropped/processed-frame counts, `dry_run`).
- **ppm-based m/z window** (`tof_half_width_for_ppm`): `--mz-ppm <ppm>` derives the
  vertical filter's TOF-index half-width from a mass tolerance at a reference m/z
  (`--mz-ppm-ref`, default the acquired-range midpoint), overriding
  `--mz-half-width`.
- **diaPASEF per-window MS/MS filtering** (`dia_window` module,
  `DiaWindowParams`, `in_window_mask`, `filter_per_window`, and
  `tdf::read_dia_windows` joining `DiaFrameMsMsInfo` + `DiaFrameMsMsWindows`):
  `--dia-per-window` runs the MS/MS filter independently within each isolation
  window's scan slice, so a mobility run is never fused across a window
  boundary (no cross-talk between the unrelated precursor m/z bands adjacent
  windows isolate). Same conditional default as the MS/MS gates above (on
  whenever MS/MS frames are filtered, and `--no-dia-per-window` reverts to
  whole-frame filtering). No-op on ddaPASEF. The same window reader backs the
  `--dia-window` and `--dia-ms1-window` gates.
- ddaPASEF **MS/MS denoising** (`msms` module, `MsmsFilterParams`,
  `combine_and_filter`): for each precursor, combine its fragment scans across the
  frames it was isolated in, denoise the combined spectrum (vertical + halo), and
  prune the individual scans to the surviving `(scan, TOF)`. Opt-in via
  `--denoise-msms` (off by default) with separate `--msms-*` knobs (library:
  `Stages::denoise_msms`). Unlike MS1 denoising, this modifies MS/MS spectra
  and identifications.
- Horizontal-halo filter (`halo` module, `HaloParams`,
  `horizontal_halo_keep_mask`): removes the weak m/z halo flanking bright ions
  (left/right only) by comparing each peak to the max of its surrounding box
  excluding its own TOF column, in integer `(scan, TOF index)` space. Runs after
  the vertical filter, **on by default**. Disable with `--no-halo` (library:
  `Stages::halo`).
- Public library API with two tiers: high-level `denoise` /
  `denoise_with_progress`, and low-level `FlatFrame`, `filter_once`,
  `filter_iterated`, `running_average`, and the type-2 `codec` module.
- Typed error types `DnoiseError` and `DecodeError` replace `anyhow` in the
  library. `anyhow` is now used only by the CLI binary.
- `denoise_with_progress` with a `Progress` callback, so the library no longer
  depends on any terminal-UI crate.
- `cli` cargo feature (enabled by default) gating the binary and its
  `clap`/`indicatif`/`serde`/`toml`/`anyhow` dependencies. Build with
  `--no-default-features` for a library-only dependency tree.
- **Experimental** `--frame-half-width` / `frame_half_width`: pre-average each
  MS1 frame over its `2r+1` MS1-frame neighborhood before filtering. See the
  README for its current limitations.

### Changed
- The SQLite `tdf` plumbing is now crate-private. The type-2 codec moved to the
  public `dnoise::codec` module (`dnoise::tdf::encode::*` → `dnoise::codec::*`).

[Unreleased]: https://github.com/pgarrett-scripps/dnoise/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.3.0
[0.1.0]: https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.1.0
