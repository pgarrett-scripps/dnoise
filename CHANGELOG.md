# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/pgarrett-scripps/dnoise/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.1.0
