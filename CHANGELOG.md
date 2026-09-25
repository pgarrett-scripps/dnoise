# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.0] - 2026-09-25

**Default output changes on ddaPASEF and diaPASEF.** The halo filter's default
`peak_fraction` is 0.10 (0.4.0: 0.15), which keeps more of each isotope
envelope; this changes MS1 output and, with `--denoise-msms`, MS/MS output.
Pass `--halo-peak-fraction 0.15` (or `halo_peak_fraction = 0.15` in the
config) for 0.4.0's value. Three further changes to the MS1 stages; the first two widen what is
kept near a gate edge:

- A feature kept by feature-level gating is now kept over its whole extent.
  0.4.0 cut it at 0.1 1/K0 beyond the mobility range of its inside points.
- Both MS1 gates are padded by default: `mz_pad = 3.0` Th and `im_pad = 0.015`
  1/K0 on each side (0.4.0: 0 and 0). To get 0.4.0's literal gate geometry,
  pass `--ms1-polygon-mz-pad 0 --ms1-polygon-im-pad 0 --dia-ms1-mz-pad 0
  --dia-ms1-im-pad 0` (the 0.1 1/K0 reach cannot be restored).
- The MS1 stage order is now streak filter, gate, halo (0.4.0: streak, halo,
  gate). Points the gate drops no longer act as halo references, and
  feature-level gating now links through points the halo filter removes
  afterwards. Temporal neighbor support keeps the old order.

### Added
- **New crate `dnoise-core`**: the in-memory denoising stages, split out of
  `dnoise` so other programs can embed them. It does no I/O and depends only on
  `serde` and `thiserror`: no `timsrust`, SQLite, zstd or rayon, and no crate
  with native code or a `links` key (checked in CI). The motivating user is
  sage-plus, which reads frames with its own `timsrust` 0.6.6 and `rusqlite`
  0.35 and calls the core per frame, in memory. Its entry point is
  `dnoise_core::Denoiser`: built once per run from plain frame metadata and
  the parameters, then one call per MS1 or MS/MS frame. The `dnoise` CLI and
  library call the same `Denoiser`, so there is one code path; output is
  byte-identical to before the split. `dnoise-core/README.md` documents which
  state crosses frames. `dnoise` keeps all I/O (the timsrust adapter, SQLite,
  the codec, the writer, provenance, batch, CLI) and re-exports the moved
  modules and types at their old paths (`dnoise::filter`, `dnoise::params`,
  `dnoise::FlatFrame`, ...).
- Polygon gate: with `im_pad > 0` the padded gate takes the polygon's exact m/z
  extent over the whole `±im_pad` band at every scan, instead of sampling only
  three scan lines (which dropped a thin spike or sharp vertex between them).
  After building, every scan's unpadded polygon interval is checked to lie in
  the padded gate, and the run fails with a clear message if not. Property
  tests (`tests/polygon_props.rs`) and `examples/polygon_check.rs` cover it.
- A containment test for padded diaPASEF MS1 windows: every point of a window
  and every point within the pads of it is kept.
- A warning when `analysis.tdf` has more than one `GroupProperties` row for the
  selection polygon; the gate uses the first.
- Each run logs (info) the raw-unit parameters of its active stages (TOF
  indices, scans) with their physical equivalents for that run's calibration,
  one line each: TOF half-widths in ppm at m/z 400, 800 and 1200, scan counts in
  1/K0 at the mid scan. The same values go into `dnoise.provenance.json` as
  `unit_equivalents` (`dnoise::units`). The parameters themselves stay in raw
  units.

### Changed
- Halo filter: default `peak_fraction` 0.10 (was 0.15). A box-kernel sweep on
  the benchmark data put 0.10 within about 1,000 identifications of a trained
  kernel at matched removal.
- MS1 m/z pre-cut: with one MS1 gate active, points far outside the gate's TOF
  range skip the streak filter. Output is unchanged (the margin covers the
  streak filter's reach, and a feature that escapes the range reruns the frame
  whole); it saves the filter work on the m/z range the gate drops anyway.
- `ms1_polygon_mz_pad` / `dia_ms1_mz_pad` default to 3.0 Th and
  `ms1_polygon_im_pad` / `dia_ms1_im_pad` to 0.015 1/K0 (were 0). MS1 gates only.
  The m/z pads are documented in Th (m/z units), not Da.
- timsrust 0.4.2 -> `timsrust-tdf` 0.6.6 behind a local adapter (`dnoise::tsr`)
  that keeps 0.4.2's converters, metadata parsing and frame order;
  `rusqlite` 0.32 -> 0.35. See `PORT_NOTES.md`.

### Deprecated
- `ms1_polygon_overlap_reach` / `dia_ms1_overlap_reach` and
  `--ms1-polygon-overlap-reach` / `--dia-ms1-overlap-reach`: accepted so 0.4.0
  configs still load, ignored with a warning, never written to the recipe.

### Removed
- `dnoise::FrameCtx` and `dnoise::process_frame_decoded`. Neither could be used
  outside the crate (`FrameCtx` had private fields and no constructor); use
  `RunContext::process` or `dnoise_core::Denoiser`. `dnoise::DecodedFrame` is
  now `dnoise_core::DecodedFrame`, re-exported at the same path.
- `overlap_reach` from `Ms1PolygonParams` and `DiaMs1WindowParams`, and the
  `reach` argument of `overlap::extend_to_features`. `PolygonGate::overlap` and
  `DiaMs1Gate::overlap` are now `bool`.


## [0.4.0] - 2026-09-23

The MS1 acquisition gates decide per feature instead of per point, and their
fixed pads are gone.

**Default MS1 output changes on ddaPASEF and diaPASEF.** 0.3.0 kept an MS1 point
only if it fell inside the selection polygon (ddaPASEF) or an isolation window
(diaPASEF), each widened by 5 Da and 0.05 1/K0. The pads were there so an edge
precursor kept its isotopes and mobility spread, but they also let back in
everything else lying near the edge, including a strip of the singly charged
band. 0.4.0 groups the points that survived the streak and m/z-halo filters into
features, using the streak filter's own adjacency, and keeps a feature whole when
any of its points lies inside the gate. A kept feature may extend at most 0.1 1/K0
beyond the mobility range of its inside points, so a long constant-m/z line
that clips the gate is not carried across the whole mobility range; bright
precursor features in the benchmark runs span under 0.1 1/K0. On a 5-minute
ddaPASEF benchmark run the most intense MS1 frame keeps 37.6% of its points
instead of 40.3%, and the MS1 area of the precursors the instrument selected
for fragmentation is retained exactly as before.

**The MS1 gates convert mobility on Bruker's calibrated scale.** Earlier
versions turned scans into 1/K0 with timsrust's straight line between the
acquisition-range bounds, but the selection polygon and the isolation windows
are defined on the run's `TimsCalibration` scale, which differs by up to
0.03 1/K0. The gates, their 1/K0 pads and reach, and the mobility crop now use
that calibration (ported from koth and checked against the timsdata SDK). On the
most intense MS1 frame of the 5-minute ddaPASEF run, 652 of 355,179 points
(0.18%) change.

The diaPASEF MS1 window gate now covers exactly a window's
`[ScanNumBegin, ScanNumEnd)` scans when `dia_ms1_im_pad` is 0. Before, the end
scan was included and float noise in the scan -> 1/K0 -> scan round trip could
add a scan at either edge. `Calibration::scan_to_im` falls back to the linear
scale with a warning when a run's calibration cannot be read; any gate that
applies to the run still refuses it.

The MS/MS isolation-window gate (`--dia-window`) is unchanged: outside a
window's scans the quadrupole passes none of its precursors.

### Added
- `ms1_polygon_overlap` / `dia_ms1_overlap` (on) and `--no-ms1-polygon-overlap`
  / `--no-dia-ms1-overlap` to gate point by point.
- `ms1_polygon_overlap_reach` / `dia_ms1_overlap_reach` (1/K0, default 0.1;
  0 = unlimited).
- `mobility_scale` (`"calibrated"`, default, or `"linear"`) and
  `--linear-mobility`; `dnoise::mobility` with the `TimsCalibration` ModelType 2
  model. `Calibration::scan_to_im` follows the run's scale.

### Changed
- `ms1_polygon_mz_pad`, `ms1_polygon_im_pad`, `dia_ms1_mz_pad` and
  `dia_ms1_im_pad` default to 0 (were 5 Da and 0.05 1/K0). To reproduce 0.3.0,
  pass `--linear-mobility --no-ms1-polygon-overlap --ms1-polygon-mz-pad 5
  --ms1-polygon-im-pad 0.05` (and the `dia-ms1` equivalents). The diaPASEF
  window gate then matches 0.3.0 to within one scan at each window's high-scan
  edge, because the pad is now measured from the window's last scan.

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

[Unreleased]: https://github.com/pgarrett-scripps/dnoise/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.5.0
[0.4.0]: https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.4.0
[0.3.0]: https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.3.0
[0.1.0]: https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.1.0
