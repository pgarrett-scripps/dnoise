# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **ddaPASEF MS/MS out-of-window gate** (`--dda-window`, `dda_window_scan_pad`,
  `DdaWindowParams`, `tdf::DiaWindows::from_pasef`): the ddaPASEF twin of
  `--dia-window`, dropping MS/MS points whose mobility scan falls outside every
  `PasefFrameMsMsInfo` isolation event for their frame. Standard timsTOF
  ddaPASEF acquisitions record MS/MS scans only inside scheduled isolation
  events, so on such files the gate removes nothing — it enforces the invariant
  rather than trusting the acquisition. Same conditional default as
  `--dia-window`: on whenever MS/MS frames are filtered, `--no-dda-window`
  forces it off. No effect on diaPASEF. `Stages` gains
  `dda_window: Option<&DdaWindowParams>`.
- **Calibration-segment warning** (`tdf::count_calibration_segments`): when an
  acquisition gate is active (`ms1_polygon` / `dia_ms1_window`) and the run's
  frames reference more than one `MzCalibration` or `TimsCalibration` segment,
  dnoise now logs a warning (in both the file writer and `RunContext`). The
  gates convert their physical-unit definitions to index space once with the
  run-level calibration timsrust exposes, which is exact only for a
  single-segment run; on a multi-segment file the gate boundary would be
  ppm-scale offset on frames referencing other segments, and the warning says
  so instead of gating silently.
- **Pinned toolchain** (`rust-toolchain.toml`, 1.97.1): the exact toolchain the
  published paper binary was built with; `rust-version = "1.85"` remains the
  MSRV.
- **Codec edge-case tests** (cargo-mutants driven): exact-bytes coverage for
  `encode_empty_frame_type2`, a round-trip with values spanning all four byte
  planes, and rejection of malformed records (short / over- / under-declared
  `total_byte_count`).

### Changed
- **Acquisition-aware noise gates are now on by default, replacing `--preset`.**
  The MS1 selection-polygon gate (`--ms1-polygon`, ddaPASEF/PASEF) and the
  diaPASEF MS1 out-of-window gate (`--dia-ms1-window`) default on; the MS/MS
  gates (`--dia-window`, `--dda-window`, `--dia-per-window`) default on whenever
  MS/MS frames are filtered (`--denoise-msms` / `--all-frames`). Each gate
  auto-detects its defining geometry and is a silent no-op when it is absent,
  so a single default set picks the right gate per acquisition. New `--no-*`
  flags (`--no-ms1-polygon`, `--no-dia-ms1-window`, `--no-dia-window`,
  `--no-dda-window`, `--no-dia-per-window`) or `<key> = false` in the config
  force a gate off. The `--preset none|auto|dda|dia` bundle (whose `auto` mode
  approximated exactly this behavior) is removed; `--preset none` callers
  should pass the `--no-*` flags instead.
- **`Ms1PolygonParams` pad defaults now match the diaPASEF MS1 gate:**
  `ms1_polygon_mz_pad` 0.0 → 5.0 Da and `ms1_polygon_im_pad` 0.0 → 0.05 1/K0,
  the same edge leniency `dia_ms1_mz_pad`/`dia_ms1_im_pad` have always applied.
  The two MS1 gates answer the same question ("could this point ever have been
  selected as a precursor?"), so they now share the same answer at the polygon
  edge: an edge precursor keeps its isotopic envelope (m/z runs high) and its
  mobility spread instead of being clipped by the literal polygon. Set both
  pads to `0` to reproduce the old behavior.
- **MS/MS streak-filter defaults retuned:** `msms_min_feature_length` 2 → 3 and
  `msms_max_internal_gap` 5 → 8, the configuration benchmarked in the paper.

### Fixed
- The example `dnoise.toml` had drifted from the code defaults
  (`max_internal_gap` was shown as 1; the built-in default is 2).

## [0.1.0] - 2026-08-03

First public release.

### Added
- **Region-of-interest crop / trim** (`crop` module, `CropParams`, `CropGate`):
  keep only points inside an axis-aligned `(m/z, 1/K0, intensity)` box plus a
  retention-time window, to carve a smaller `.d` out of a large one. CLI:
  `--mz-min`/`--mz-max` (Da), `--im-min`/`--im-max` (1/K0),
  `--min-intensity`/`--max-intensity`, and `--rt-min`/`--rt-max` (minutes). m/z and
  mobility become integer `(TOF, scan)` ranges via the run calibration and apply to
  every frame; RT bounds empty out-of-window frames (never deleted, so the frame
  axis and dependent tables stay valid and SDK-compatible). Composes with the
  denoiser, or run `--crop-only` to apply just the crop with no denoising.
- **Named presets & acquisition auto-detection** (`Acquisition`,
  `detect_acquisition`): `--preset auto` reads the run's `MsMsType` and enables the
  MS1 gate(s) matched to the scheme (ddaPASEF → selection-polygon; diaPASEF →
  isolation-window gates); `--preset dda`/`--preset dia` force a bundle. Presets
  supply gate defaults only; explicit gate flags and config values still win.
- **Dry runs, sampling & JSON reports**: `--dry-run` runs the full pipeline and
  reports the reduction without writing any output; `--sample <f>` (dry-run only)
  processes a deterministic fraction of frames for a fast estimate
  (`--sample-seed`); `--report <FILE>` writes the effective config plus per-MS-level
  reduction statistics as JSON. Library: `denoise_with_options`, `RunOptions`,
  `SampleSpec`, and a widened `DenoiseStats` (per-level point counts, summed
  intensities, cropped/processed-frame counts, `dry_run`).
- **ppm-based m/z window** (`tof_half_width_for_ppm`): `--mz-ppm <ppm>` derives the
  vertical filter's TOF-index half-width from a mass tolerance at a reference m/z
  (`--mz-ppm-ref`, default the acquired-range midpoint), overriding
  `--mz-half-width`.
- **diaPASEF isolation-window MS/MS filtering** (`dia_window` module,
  `DiaWindowParams`, `in_window_mask`, `filter_per_window`; `tdf::read_dia_windows`
  joining `DiaFrameMsMsInfo` + `DiaFrameMsMsWindows`). Two opt-in, diaPASEF-only
  features sharing one window reader: `--dia-window` drops MS/MS points whose
  mobility scan falls outside every isolation window for their frame
  (out-of-window noise; `--dia-window-scan-pad` adds edge leniency), and
  `--dia-per-window` runs the MS/MS filter independently within each window's scan
  slice so a mobility run is never fused across a window boundary (no cross-talk
  between the unrelated precursor m/z bands adjacent windows isolate). Both are
  no-ops on ddaPASEF. `denoise`/`denoise_with_progress` gain
  `dia_window: Option<&DiaWindowParams>` and `dia_per_window: bool` arguments.
- ddaPASEF **MS/MS denoising** (`msms` module, `MsmsFilterParams`,
  `combine_and_filter`): for each precursor, combine its fragment scans across the
  frames it was isolated in, denoise the combined spectrum (vertical + halo), and
  prune the individual scans to the surviving `(scan, TOF)`. Opt-in via
  `--denoise-msms` (off by default) with separate `--msms-*` knobs;
  `denoise`/`denoise_with_progress` gain a `denoise_msms: Option<&MsmsFilterParams>`
  argument. Unlike MS1 denoising, this modifies MS/MS spectra and identifications.
- Horizontal-halo filter (`halo` module, `HaloParams`,
  `horizontal_halo_keep_mask`): removes the weak m/z halo flanking bright ions
  (left/right only) by comparing each peak to the max of its surrounding box
  excluding its own TOF column, in integer `(scan, TOF index)` space. Runs after
  the vertical filter, **on by default**; disable with `--no-halo`.
  `denoise`/`denoise_with_progress` gain a
  `halo: Option<&HaloParams>` argument.
- Public library API with two tiers: high-level `denoise` /
  `denoise_with_progress`, and low-level `FlatFrame`, `filter_once`,
  `filter_iterated`, `running_average`, and the type-2 `codec` module.
- Typed error types `DnoiseError` and `DecodeError` (replacing `anyhow` in the
  library; `anyhow` is now used only by the CLI binary).
- `denoise_with_progress` with a `Progress` callback, so the library no longer
  depends on any terminal-UI crate.
- `cli` cargo feature (enabled by default) gating the binary and its
  `clap`/`indicatif`/`serde`/`toml`/`anyhow` dependencies. Build with
  `--no-default-features` for a library-only dependency tree.
- **Experimental** `--frame-half-width` / `frame_half_width`: pre-average each
  MS1 frame over its `2r+1` MS1-frame neighborhood before filtering. See the
  README for its current limitations.

### Changed
- The SQLite `tdf` plumbing is now crate-private; the type-2 codec moved to the
  public `dnoise::codec` module (`dnoise::tdf::encode::*` → `dnoise::codec::*`).

[Unreleased]: https://github.com/pgarrett-scripps/dnoise/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/pgarrett-scripps/dnoise/releases/tag/v0.1.0
