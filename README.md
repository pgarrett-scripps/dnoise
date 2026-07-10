# dnoise

[![Crates.io](https://img.shields.io/crates/v/dnoise.svg)](https://crates.io/crates/dnoise)
[![Docs.rs](https://docs.rs/dnoise/badge.svg)](https://docs.rs/dnoise)
[![License](https://img.shields.io/crates/l/dnoise.svg)](#license)

Denoise Bruker timsTOF `.d` folders by reducing the raw 3-D data
(`scan × TOF-index × intensity`) with the **iterative vertical-IM feature
filter** (Stage 1 of [ALGORITHM.md](ALGORITHM.md)). Centroiding and smoothing
are intentionally not implemented.

A real ion forms a vertical streak in `(TOF-index × scan)` space. The filter
walks each TOF index, sums the ion-mobility profile in a small TOF window
around it, and keeps only points belonging to long-enough, intense-enough
vertical runs. The pass is iterated: each pass re-filters the survivors of the
previous one.

`dnoise` is both a **command-line tool** and a **Rust library**.

## Install

The CLI:

```bash
cargo install dnoise
```

As a library dependency (no CLI deps pulled in):

```toml
[dependencies]
dnoise = { version = "0.1", default-features = false }
```

The `cli` feature (on by default) adds the binary and its `clap`/`indicatif`/
`serde`/`toml`/`anyhow` dependencies; disable default features for library-only
use.

## Command-line usage

```bash
dnoise <INPUT.d> <OUTPUT.d> [options]
dnoise <INPUT.d> --in-place [options]
```

By default the source folder is never modified; a new `.d` is written with a rewritten
`analysis.tdf_bin` (re-encoded as compression **type 2**) and an updated
`analysis.tdf` (`Frames.TimsId/NumPeaks/MaxIntensity/SummedIntensities` and
`GlobalMetadata.TimsCompressionType`). The leading reserved header that Bruker
places before the first frame (the smallest `Frames.TimsId`, often 64 bytes) is
copied verbatim and all rewritten offsets are shifted past it, so the output is
byte-layout-compatible with the Bruker SDK / `timsdata` DLL.

| Option | Default | Meaning |
|---|---|---|
| `--mz-half-width` | 3 | Column half-width in TOF indices (`[c-w, c+w]`). |
| `--min-feature-length` | 5 | Minimum total span (scans) of a kept feature. |
| `--max-internal-gap` | 2 | Max empty scans tolerated inside a feature. |
| `--min-window-intensity` | 0 | Per-scan summed-intensity floor for occupancy. |
| `--min-feature-intensity` | 0 | Total summed-intensity floor for a kept feature. |
| `--iterations` | 2 | Filter passes (each re-applies to prior survivors). |
| `--no-halo` | (on) | Disable the horizontal-halo filter, which runs after the vertical filter (see below). |
| `--halo-peak-fraction` | 0.1 | Drop a peak below this fraction of its off-column box-max. |
| `--halo-mz-idx-half-width` | 100 | Reference-box half-width along TOF index. |
| `--halo-scan-half-width` | 2 | Reference-box half-width along ion-mobility scan. |
| `--denoise-msms` | off | Denoise ddaPASEF **MS/MS** frames precursor-by-precursor (see below). Changes MS/MS spectra and IDs. |
| `--msms-min-feature-length` | 3 | MS/MS filter: min vertical-run span (+ `--msms-mz-half-width`, `--msms-max-internal-gap`, `--msms-min-window-intensity`, `--msms-min-feature-intensity`, `--msms-iterations`). |
| `--dia-window` | off | **diaPASEF only.** Drop MS/MS points whose mobility scan falls outside every isolation window for their frame (see below). |
| `--dia-window-scan-pad` | 0 | Scans of leniency added to each side of every isolation window before a point counts as out-of-window. |
| `--dia-per-window` | off | **diaPASEF only.** When the MS/MS filter runs, filter each isolation window's scan slice independently (no cross-window linking; see below). |
| `--dia-ms1-window` | off | **diaPASEF only.** Drop MS1 points whose `(m/z, mobility)` falls outside every isolation window (precursors that are never fragmented; see below). |
| `--dia-ms1-mz-pad` | 5 | MS1 gate: m/z leniency (Da) added to each side of every window, so an edge precursor keeps its full isotopic envelope. |
| `--dia-ms1-im-pad` | 0.05 | MS1 gate: ion-mobility leniency (1/K0) added to each side of every window. |
| `--frame-half-width` | 0 | **Experimental.** Pre-average each MS1 frame over its `2r+1` MS1-frame neighborhood before filtering (see below). |
| `--all-frames` | off | Also filter MS/MS frames (default: MS1 only). |
| `--preset` | `none` | Bundle of MS1 gates for the acquisition scheme: `auto` detects ddaPASEF/diaPASEF and enables the matching gate(s), `dda` = selection-polygon, `dia` = isolation-window gates, `none` = nothing. Explicit gate flags still win (see below). |
| `--mz-ppm` | — | Set the vertical filter's m/z window from a mass tolerance in ppm (converted at a reference m/z) instead of raw TOF indices; overrides `--mz-half-width` (see below). |
| `--mz-ppm-ref` | acq midpoint | Reference m/z (Da) for `--mz-ppm`. |
| `--mz-min` / `--mz-max` | — | **Crop.** Keep only points in this m/z band (Da). |
| `--im-min` / `--im-max` | — | **Crop.** Keep only points in this ion-mobility band (1/K0). |
| `--rt-min` / `--rt-max` | — | **Crop.** Keep only frames in this retention-time window (**minutes**); out-of-window frames are emitted empty (see below). |
| `--min-intensity` / `--max-intensity` | — | **Crop.** Keep only points in this intensity range. |
| `--crop-only` | off | Apply only the crop and skip all denoising, so the output is a raw subset of the input (see below). |
| `--dry-run` | off | Estimate the reduction without writing any output (see below). |
| `--sample` | — | With `--dry-run`, process only this fraction (`0 < f ≤ 1`) of frames, chosen deterministically, for a fast estimate. |
| `--sample-seed` | 0 | Seed for `--sample` frame selection. |
| `--report` | — | Write a JSON run report (effective config + reduction stats) to this file. |
| `--threads` | all cores | Worker threads. |
| `--config` / `-c` | — | Load parameters from a TOML file (see below). |
| `--force` | off | Overwrite an existing output folder. |
| `--in-place` | off | Denoise the input folder in place (omit `OUTPUT`); see below. |

Filtering runs in parallel across frames (rayon); frames are written in order
so binary offsets stay consistent.

> **Presets & auto-detection (`--preset`).** The right MS1 gate depends on how the
> run was acquired: ddaPASEF benefits from the selection-polygon gate, diaPASEF from
> the isolation-window gates. `--preset auto` reads the run's `MsMsType` and enables
> the matching gate(s) for you (ddaPASEF → `--ms1-polygon`; diaPASEF →
> `--dia-ms1-window --dia-window`); MS1-only or unrecognised schemes enable nothing.
> `--preset dda` / `--preset dia` force a bundle regardless of detection. A preset
> only supplies gate *defaults*: an explicit gate flag (or config value) always wins,
> and gates that find nothing to act on (e.g. a run with no polygon) are silently
> skipped as usual.

> **Region-of-interest crop (`--mz-*` / `--im-*` / `--rt-*` / `--*-intensity`).** A
> crop is a blunt subset of the raw acquisition (not a signal/noise decision): it is
> how you carve a smaller `.d` out of a large one for sharing, faster downstream
> searches, or test fixtures. m/z and mobility bounds become integer `(TOF, scan)`
> ranges via the run calibration and apply to **every** frame (MS1 and MS/MS alike);
> the intensity range is a per-point floor/ceiling. Retention-time bounds act at the
> frame level: an out-of-window frame is emitted **empty** rather than deleted, so
> the frame axis and every table that references it stay valid and Bruker-SDK
> compatible. The crop composes with the denoiser (applied as an extra filter), or
> run it alone with `--crop-only` to leave retained signal untouched. Note the crop
> does not rewrite the acquisition-range metadata (`MzAcqRange…`), which continues to
> describe what the instrument acquired.

> **Dry runs & reports (`--dry-run`, `--sample`, `--report`).** `--dry-run` runs the
> full pipeline but writes nothing, printing the reduction so you can tune parameters
> without producing an output `.d`. Add `--sample 0.1` to process a deterministic 10%
> of frames for a fast estimate (the ratio is representative; only valid with
> `--dry-run`). `--report out.json` writes the effective configuration plus the
> reduction statistics (per-MS-level point counts, summed intensities, cropped-frame
> count, elapsed time) as JSON, for parameter sweeps or provenance; it works for both
> real and dry runs.

> **ppm-based m/z window (`--mz-ppm`).** `--mz-half-width` is a constant in TOF
> indices, but a real peak's width scales with m/z. `--mz-ppm 20` derives the window
> from a mass tolerance instead, evaluated at a reference m/z (`--mz-ppm-ref`, default
> the midpoint of the acquired range) via the run calibration. The vertical filter
> still uses one constant index window, so this sets a physically meaningful width
> once rather than making the window vary across the mass range.

> **In-place denoising (`--in-place`).** Omit the `OUTPUT` argument and pass
> `--in-place` to overwrite the input folder. dnoise still never edits the source
> while reading it: it writes a temporary sibling folder (`<INPUT>.dnoise-tmp`)
> and, only after a clean run, moves the original aside to `<INPUT>.dnoise-old`,
> renames the new folder into place, and deletes the backup. If the final rename
> fails the original is restored, so an interrupted run leaves the input intact
> (plus a recoverable `.dnoise-tmp`/`.dnoise-old` sibling). The two renames are
> atomic only when the temp folder lands on the same filesystem as the input,
> which it does by construction.

> **MS/MS denoising (ddaPASEF, opt-in via `--denoise-msms`).** ddaPASEF
> fragments a precursor across one ion-mobility scan window repeated over several
> frames. dnoise combines each precursor's fragment scans across those frames into
> one spectrum (summing intensity at aligned `(scan, TOF)`), runs the vertical +
> halo filters on the combined spectrum, and then prunes the individual scans to
> the surviving `(scan, TOF)`. Unlike MS1 denoising this **modifies MS/MS spectra
> and therefore identifications** — measure its effect by re-searching. Tuned with
> the `--msms-*` knobs (separate from the MS1 ones; default `min_feature_length`
> is smaller because the windows are short).

> **diaPASEF isolation windows (opt-in, diaPASEF only).** In diaPASEF the
> quadrupole steps through a set of `(mobility, m/z)` isolation windows per cycle;
> each window occupies a contiguous mobility-scan interval (`DiaFrameMsMsWindows`).
> Two opt-in features use that scheme. `--dia-window` drops MS/MS points whose scan
> falls outside every window for their frame — signal that was never isolated
> (out-of-window noise, typically at the mobility edges); `--dia-window-scan-pad`
> widens each window to tolerate signal just past an edge. `--dia-per-window` makes
> the MS/MS filter (`--denoise-msms` or `--all-frames`) run **independently inside
> each window's scan slice** instead of over the whole frame, so the vertical
> filter cannot fuse a mobility run across a window boundary — i.e. no cross-talk
> between the unrelated precursor m/z bands that adjacent windows isolate. Both are
> no-ops on ddaPASEF (no `DiaFrameMsMs*` tables).

> **diaPASEF MS1 out-of-window gate (`--dia-ms1-window`, opt-in, diaPASEF only).**
> The union of all isolation windows is the precursor space the method can ever
> fragment. This gate drops **MS1** points whose `(m/z, mobility)` falls in no
> window — precursors that are never selected — keeping the survey scans to the
> useful precursor band. Each window is padded in **physical units**:
> `--dia-ms1-mz-pad` (Da, default 5) and `--dia-ms1-im-pad` (1/K0, default 0.05),
> converted to TOF indices / scans once via the run's calibration, so a precursor
> sitting at a window edge keeps its full isotopic envelope (isotopes run to higher
> m/z) and mobility spread. Applies to MS1 frames only and is a no-op on ddaPASEF.

> **Horizontal-halo filter (on by default).** After the vertical filter, dnoise
> removes the weak m/z halo flanking bright ions — left/right only. Each peak is
> compared to the maximum intensity in its surrounding box (`±halo-scan-half-width`
> scans × `±halo-mz-idx-half-width` TOF indices) **excluding its own TOF column**,
> and dropped if its intensity is below `peak_fraction` of that reference.
> Excluding the own column means the vertical streak above/below a peak never
> counts against it — only genuine left/right neighbors do. It works in integer
> `(scan, TOF index)` space (no calibration) and keeps/drops native points (no
> smoothing). Ported from `tdfpy`'s `HorizontalHaloFilter`; disable with
> `--no-halo`.

> **`--frame-half-width` is experimental.** It replaces each MS1 frame with the
> centered running average of its `2r+1` MS1-frame neighborhood before
> filtering. With the default zero intensity thresholds the exact-`(scan, tof)`
> merge mostly *concatenates* adjacent frames (they share only ~4% of bins) and
> inflates the output rather than denoising it; it only suppresses noise when
> `--min-window-intensity` is raised above the single-frame noise floor. Leave
> it at `0` unless you are deliberately experimenting.

### Config file

Instead of (or alongside) flags, parameters can come from a TOML file:

```bash
dnoise <INPUT.d> <OUTPUT.d> --config dnoise.toml
```

Every key is optional and uses the same name as the flag with underscores. See
[dnoise.toml](dnoise.toml) for a fully-commented example. Precedence is
**explicit CLI flag > config file > built-in default**, so a config sets the
baseline and individual flags override it for one run. Unknown keys are rejected
to catch typos.

```toml
mz_half_width = 3
min_feature_length = 7
iterations = 2
all_frames = false
# threads = 8
```

### Logging

`dnoise` narrates each run on **stderr** and prints only the final result line on
**stdout**, so the two streams can be captured independently (handy for scripts and
AI tooling):

```bash
dnoise <INPUT.d> <OUTPUT.d> 2> run.log      # logs to run.log, result to the console
```

The logs are structured, one event per line, and cover the effective configuration
(every resolved knob and which stages are enabled), the detected acquisition scheme
(ddaPASEF / diaPASEF / MS1-only) and frame inventory, which gates activated or were
skipped and why, progress, and a final `denoise: complete` with raw/kept point
counts. Example (abridged):

```text
INFO dnoise: config: enabled stages halo=true ms1_polygon=true denoise_msms=false ...
INFO dnoise::writer: denoise: frame inventory scheme="ddaPASEF" frames=8639 ms1=786 msms=7853 empty=0
INFO dnoise::writer: MS1 selection-polygon gate active
INFO dnoise::writer: denoise: complete frames=8639 raw_points=300509979 kept_points=110092035 kept_pct=36.64
```

Verbosity: `-v` adds debug detail (e.g. why a requested gate was skipped), `-vv`
adds trace; `-q` limits output to warnings and errors. For fine-grained control set
`RUST_LOG` (e.g. `RUST_LOG=dnoise=debug`), which overrides the flags. When stderr is
an interactive terminal a progress bar is shown instead of periodic progress lines.

## Library usage

Two API tiers are exposed (full docs on [docs.rs](https://docs.rs/dnoise)):

**High level** — denoise a whole `.d` folder:

```rust,no_run
use dnoise::{FilterParams, denoise};
use std::path::Path;

let stats = denoise(
    Path::new("input.d"),
    Path::new("output.d"),
    &FilterParams::default(),
    false, // filter MS1 only
    0,     // frame_half_width (0 = no pre-averaging)
    false, // don't overwrite an existing output
)?;
println!("{} -> {} points", stats.raw_points, stats.kept_points);
# Ok::<(), dnoise::DnoiseError>(())
```

Use `denoise_with_progress` to receive `Progress` updates (the CLI uses this to
drive its progress bar) — the library itself depends on no UI crate.

**Low level** — run the filter on an in-memory frame, or use the type-2
[`codec`](https://docs.rs/dnoise/latest/dnoise/codec/) directly:

```rust
use dnoise::{FilterParams, FlatFrame, filter_iterated};

let frame = FlatFrame {
    frame_id: 1,
    num_scans: 700,
    scan: (10..18).collect(),   // a streak across 8 mobility scans
    tof: vec![1000; 8],
    intensity: vec![100; 8],
};
let keep = filter_iterated(&frame, &FilterParams::default());
assert!(keep.iter().all(|&k| k));
```

## Validate output

```bash
cargo run --release --example validate -- <PATH.d>
```

Re-reads every frame with timsrust and checks each frame's decoded peak count
against `Frames.NumPeaks`.

To verify the type-2 codec against *real Bruker bytes* (decode raw frames
straight from `analysis.tdf_bin` and compare to timsrust):

```bash
cargo run --release --example check_codec -- <PATH.d> [num_frames]
```

## Limitation: input compression type

dnoise reads **compression type 2** (and uncompressed) input and always writes
type 2. It cannot yet read **type 3** (zstd + bitshuffle) `.d` files: timsrust's
`timscompress` feature depends on a `timscompress` crate that is only a stub on
crates.io, so the decoder does not build. Type-3 support is blocked on that
crate being published.

## Tests

```bash
just check     # clippy -D warnings + fmt check + tests
```

`tests/roundtrip.rs` gates the type-2 encoder by decoding what it encodes;
`tests/filter.rs` checks streak retention, gap-closing, window aggregation, and
iteration monotonicity; `tests/average.rs` covers the running-average pre-pass.

## Reproducing the paper

The benchmark suite, manuscript, and Supporting Information that accompany the
dnoise paper live under [`benchmark/`](benchmark/) and [`paper/`](paper/). The
exact state used for the paper is frozen on the [`paper`](../../tree/paper)
branch, so every figure, table, and PDF can be regenerated as published while
`main` tracks ongoing development. Start from
[`benchmark/README.md`](benchmark/README.md), which lists, in order, the commands
that rebuild each manuscript asset from the raw `.d` files (PRIDE
[PXD070049](https://www.ebi.ac.uk/pride/archive/projects/PXD070049)). Rust deps
are pinned in `Cargo.lock`; the Python analysis env is pinned in
`benchmark/pyproject.toml` + `benchmark/uv.lock` (run via `uv`).

## License

Licensed under the [MIT License](LICENSE-MIT).
