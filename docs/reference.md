# dnoise reference

Complete documentation of every CLI option, gate, and stage. For a quick
start see the [README](../README.md), for the algorithm itself see
[ALGORITHM.md](../ALGORITHM.md), and for the library API see
[docs.rs](https://docs.rs/dnoise).

For prm-PASEF, the ordinary command preserves fragments; `--denoise-msms` enables
experimental filtering within each PRM isolation event using the `--msms-*`
parameters. Discovery gates remain disabled. PRM rejects `--all-frames`, while
mixed/unknown acquisitions reject both fragment-filter options. Explicit cropping
still applies to MS/MS. PRM metadata checks remain active with `--skip-validation`.
See [targeted support](targeted.md).

## Command line

```bash
dnoise <INPUT.d> <OUTPUT.d> [options]
dnoise <INPUT.d> --in-place [options]
```

By default the source folder is never modified. A new `.d` is written with a rewritten
`analysis.tdf_bin` (re-encoded as compression **type 2**) and an updated
`analysis.tdf` (`Frames.TimsId/NumPeaks/MaxIntensity/SummedIntensities` and
`GlobalMetadata.TimsCompressionType`). The leading reserved header that Bruker
places before the first frame (the smallest `Frames.TimsId`, often 64 bytes) is
copied verbatim and all rewritten offsets are shifted past it, so the output is
byte-layout-compatible with the Bruker SDK / `timsdata` DLL.

Filtering runs in parallel across frames (rayon). Frames are written in order
so binary offsets stay consistent.

## All options

| Option | Default | Meaning |
|---|---|---|
| `--mz-half-width` | 3 | Column half-width in TOF indices (`[c-w, c+w]`). |
| `--min-feature-length` | 5 | Minimum total span (scans) of a kept feature. |
| `--max-internal-gap` | 2 | Max empty scans tolerated inside a feature. |
| `--min-window-intensity` | 0 | Per-scan summed-intensity floor for occupancy. |
| `--min-feature-intensity` | 0 | Total summed-intensity floor for a kept feature. |
| `--iterations` | 2 | Filter passes (each re-applies to prior survivors). |
| `--no-halo` | (on) | Disable the horizontal-halo filter, which runs after the vertical filter (see below). |
| `--halo-peak-fraction` | 0.15 | Drop a peak below this fraction of its off-column box-max. |
| `--halo-mz-idx-half-width` | 80 | Reference-box half-width along TOF index. |
| `--halo-scan-half-width` | 2 | Reference-box half-width along ion-mobility scan. |
| `--denoise-msms` | off | Denoise ddaPASEF **MS/MS** frames precursor-by-precursor (see below). Changes MS/MS spectra and IDs. |
| `--msms-min-feature-length` | 3 | MS/MS filter: min vertical-run span (+ `--msms-mz-half-width`, `--msms-max-internal-gap`, `--msms-min-window-intensity`, `--msms-min-feature-intensity`, `--msms-iterations`). |
| `--dia-window` | conditional | **diaPASEF only.** Drop MS/MS points whose mobility scan falls outside every isolation window for their frame. On whenever MS/MS frames are filtered (`--denoise-msms` / `--all-frames`). `--no-dia-window` forces it off (see below). |
| `--dia-window-scan-pad` | 0 | Scans of leniency added to each side of every isolation window before a point counts as out-of-window. |
| `--dia-per-window` | conditional | **diaPASEF only.** When the MS/MS filter runs, filter each isolation window's scan slice independently (no cross-window linking, see below). On whenever MS/MS frames are filtered. `--no-dia-per-window` reverts to whole-frame filtering. |
| `--dda-window` | conditional | **ddaPASEF only.** Drop MS/MS points whose mobility scan falls outside every `PasefFrameMsMsInfo` isolation event for their frame. Standard acquisitions record no such points, so this enforces an invariant rather than reducing data. On whenever MS/MS frames are filtered. `--no-dda-window` forces it off. |
| `--dda-window-scan-pad` | 0 | Scans of leniency added to each side of every isolation event. |
| `--dia-ms1-window` | on | **diaPASEF only.** Drop MS1 points whose `(m/z, mobility)` falls outside every isolation window (precursors that are never fragmented, see below). `--no-dia-ms1-window` disables. |
| `--dia-ms1-overlap` | on | MS1 gate: keep a whole streak-filter feature when any of its points lies in a window. `--no-dia-ms1-overlap` gates point by point (0.3.0). |
| `--dia-ms1-mz-pad` | 3 | MS1 gate: m/z leniency (Th) added to each side of every window (0.4.0: 0; 0.3.0: 5). |
| `--dia-ms1-im-pad` | 0.015 | MS1 gate: ion-mobility leniency (1/K0) added to each side of every window (0.4.0: 0; 0.3.0: 0.05). |
| `--ms1-polygon` | on | **ddaPASEF.** Drop MS1 points outside the run's PASEF selection polygon (never-selected precursor space). Auto-detected, so it is a no-op if the run stores no polygon or defines a diaPASEF window scheme. `--no-ms1-polygon` disables. Feature-level gating: `--ms1-polygon-overlap` (on; `--no-ms1-polygon-overlap` for point by point). Pads: `--ms1-polygon-mz-pad` (Th, default 3), `--ms1-polygon-im-pad` (1/K0, default 0.015). |
| `--smooth` | off | Final stage: box-average each survivor's intensity over its `(scan, TOF-index)` box (stabilises the watershed centroider). Sub: `--smooth-mz-idx-half-width`, `--smooth-scan-half-width`, `--smooth-iterations`. |
| `--watershed` | off | Final stage: watershed centroiding, collapsing point groups into intensity-weighted centroids (lossy). Sub: `--watershed-box-scan`, `--watershed-box-mz-idx`, `--watershed-min-seed-intensity`, `--watershed-min-centroid-total`, `--watershed-max-tof-offset`. |
| `--box-centroid` | off | Final stage: greedy small-box centroiding, tiling streaks into small centroids rather than collapsing them. Mutually exclusive with `--watershed`. Sub: `--box-centroid-mz-idx-half`, `--box-centroid-scan-half`, `--box-centroid-min-total`. |
| `--ms1-neighbor-radius` | 0 | **Experimental.** Previous/next compatible MS1 observations for the keep decision. Alias: `--frame-half-width`. |
| `--prm-neighbor-radius` | 0 | **Experimental.** Previous/next compatible observations of each PRM target; requires MS/MS denoising. |
| `--dia-neighbor-radius` | 0 | **Experimental.** Previous/next compatible observations of each DIA isolation window; requires MS/MS filtering. |
| `--neighbor-max-rt-gap` | 5 | Maximum absolute RT distance between a supporting observation and the current frame, in **seconds**. |
| `--all-frames` | off | Also filter MS/MS frames (default: MS1 only). |
| `--mz-ppm` | - | Set the vertical filter's m/z window from a mass tolerance in ppm (converted at a reference m/z) instead of raw TOF indices. Overrides `--mz-half-width` (see below). |
| `--mz-ppm-ref` | acq midpoint | Reference m/z (Da) for `--mz-ppm`. |
| `--mz-min` / `--mz-max` | - | **Crop.** Keep only points in this m/z band (Da). |
| `--im-min` / `--im-max` | - | **Crop.** Keep only points in this ion-mobility band (1/K0). |
| `--rt-min` / `--rt-max` | - | **Crop.** Keep only frames in this retention-time window (**minutes**). Out-of-window frames are emitted empty (see below). |
| `--min-intensity` / `--max-intensity` | - | **Crop.** Keep only points in this intensity range. |
| `--crop-only` | off | Apply only the crop and skip all denoising, so the output is a raw subset of the input (see below). |
| `--dry-run` | off | Estimate the reduction without writing any output (see below). |
| `--sample` | - | With `--dry-run`, process only this fraction (`0 < f ≤ 1`) of frames, chosen deterministically, for a fast estimate. |
| `--sample-seed` | 0 | Seed for `--sample` frame selection. |
| `--report` | - | Write a JSON run report (effective config + reduction stats) to this file. |
| `--threads` | all cores | Worker threads. |
| `--frame-batch-size` | 2048 | Frames encoded per batch; smaller values may lower memory at a runtime cost. |
| `--config` / `-c` | - | Load parameters from a TOML file (see below). |
| `--force` | off | Overwrite an existing output folder. |
| `--in-place` | off | Denoise the input folder in place (omit `OUTPUT`, see below). |

## Acquisition-aware gates (on by default)

The right gate depends on how the run was acquired, but each gate auto-detects
its defining geometry and is a silent no-op when it is absent, so a single
default set picks the right gate per acquisition. The MS1 gates
(`--ms1-polygon` for ddaPASEF/PASEF, `--dia-ms1-window` for diaPASEF) are
always on by default. The MS/MS gates (`--dia-window`, `--dda-window`,
`--dia-per-window`) default on only when MS/MS frames are actually filtered
(`--denoise-msms` or `--all-frames`), so a plain MS1-only run leaves fragment
spectra untouched. Force any gate off with its `--no-*` flag (or
`<key> = false` in the config). An explicit flag or config value always wins
over the default.

## Region-of-interest crop (`--mz-*` / `--im-*` / `--rt-*` / `--*-intensity`)

A crop is a blunt subset of the raw acquisition (not a signal/noise decision):
it is how you carve a smaller `.d` out of a large one for sharing, faster
downstream searches, or test fixtures. m/z and mobility bounds become integer
`(TOF, scan)` ranges via the run calibration and apply to **every** frame (MS1
and MS/MS alike). The intensity range is a per-point floor/ceiling.
Retention-time bounds act at the frame level: an out-of-window frame is
emitted **empty** rather than deleted, so the frame axis and every table that
references it stay valid and Bruker-SDK compatible. The crop composes with the
denoiser (applied as an extra filter), or run it alone with `--crop-only` to
leave retained signal untouched. Note the crop does not rewrite the
acquisition-range metadata (`MzAcqRange…`), which continues to describe what
the instrument acquired. Physical bounds reject multiple calibration references
in the dimension being cropped; RT/intensity-only crops remain available.
See [calibration restrictions](provenance.md#calibration-restrictions).

## Dry runs & reports (`--dry-run`, `--sample`, `--report`)

`--dry-run` runs the full pipeline but writes nothing, printing the reduction
so you can tune parameters without producing an output `.d`. Add
`--sample 0.1` to process a deterministic 10% of frames for a fast estimate
(the ratio is representative, and it is only valid with `--dry-run`).
`--report out.json`
writes the effective configuration plus the reduction statistics (per-MS-level
point counts, summed intensities, cropped-frame count, elapsed time) as JSON,
for parameter sweeps or provenance. It works for both real and dry runs.

## ppm-based m/z window (`--mz-ppm`)

`--mz-half-width` is a constant in TOF indices, but a real peak's width scales
with m/z. `--mz-ppm 20` derives the window from a mass tolerance instead,
evaluated at a reference m/z (`--mz-ppm-ref`, default the midpoint of the
acquired range) via the run calibration. The vertical filter still uses one
constant index window, so this sets a physically meaningful width once rather
than making the window vary across the mass range.

## In-place replacement

`dnoise input.d --in-place` processes and validates an owned temporary copy, then
replaces the input. Existing output is preserved on processing errors or cancel.
Installation failures restore the original or report its backup location.
See [file safety and crash recovery](provenance.md#file-safety-and-recovery).

## MS/MS denoising (ddaPASEF, opt-in via `--denoise-msms`)

ddaPASEF fragments a precursor across one ion-mobility scan window repeated
over several frames. dnoise combines each precursor's fragment scans across
those frames into one spectrum (summing intensity at aligned `(scan, TOF)`),
runs the vertical + halo filters on the combined spectrum, and then prunes the
individual scans to the surviving `(scan, TOF)`. Unlike MS1 denoising this
**modifies MS/MS spectra and therefore identifications**. Measure its effect
by re-searching. Tuned with the `--msms-*` knobs, separate from the MS1 ones,
with a smaller default `min_feature_length` because the windows are short.

## diaPASEF isolation windows

In diaPASEF the quadrupole steps through a set of `(mobility, m/z)` isolation
windows per cycle. Each window occupies a contiguous mobility-scan interval
(`DiaFrameMsMsWindows`). Two features use that scheme. `--dia-window` drops
MS/MS points whose scan falls outside every window for their frame: signal
that was never isolated (out-of-window noise, typically at the mobility
edges). `--dia-window-scan-pad` widens each window to tolerate signal just
past an edge. `--dia-per-window` makes the MS/MS filter (`--denoise-msms` or
`--all-frames`) run **independently inside each window's scan slice** instead
of over the whole frame, so the vertical filter cannot fuse a mobility run
across a window boundary, i.e. no cross-talk between the unrelated precursor
m/z bands that adjacent windows isolate. Both are no-ops on ddaPASEF (no
`DiaFrameMsMs*` tables).

## diaPASEF MS1 out-of-window gate (`--dia-ms1-window`, diaPASEF only)

The union of all isolation windows is the precursor space the method can ever
fragment. This gate drops **MS1** points whose `(m/z, mobility)` falls in no
window (precursors that are never selected), keeping the survey scans to the
useful precursor band. Applies to MS1 frames only and is a no-op on ddaPASEF.

The gate decides per **feature**, not per point (`--dia-ms1-overlap`, on by
default). The points that survived the streak and halo filters are grouped
with the streak filter's own adjacency (within `mz-half-width` TOF indices and
`max-internal-gap + 1` scans), and a feature is kept whole when any of its
points lies inside a window. A precursor whose mobility peak straddles a window
edge therefore keeps its full mobility spread, and a feature lying wholly
outside every window is dropped. A kept feature is kept over its whole
mobility extent. (0.4.0 capped that at `--dia-ms1-overlap-reach`, default 0.1
1/K0; the cap was removed in 0.5.0 and the `--dia-ms1-overlap-reach` /
`--ms1-polygon-overlap-reach` flags and config keys are accepted but ignored,
with a warning.)

Pads widen each window in **physical units** before the test:
`--dia-ms1-mz-pad` (Th, default 3) and `--dia-ms1-im-pad` (1/K0, default 0.015),
converted to TOF indices / scans once via the run's calibration; every scan
within the mobility pad gets the full padded m/z band. dnoise 0.4.0 used no
pads (`--dia-ms1-mz-pad 0 --dia-ms1-im-pad 0`). dnoise 0.3.0
gated point by point with pads of 5 Th and 0.05 1/K0; `--no-dia-ms1-overlap
--dia-ms1-mz-pad 5 --dia-ms1-im-pad 0.05` reproduces it. The ddaPASEF selection
polygon gate (`--ms1-polygon`) works the same way with its own
`--ms1-polygon-*` options.

## Horizontal-halo filter (on by default)

After the vertical filter, dnoise removes the weak m/z halo flanking bright
ions, left/right only. Each peak is compared to the maximum intensity in its
surrounding box (`±halo-scan-half-width` scans × `±halo-mz-idx-half-width` TOF
indices) **excluding its own TOF column**, and dropped if its intensity is
below `peak_fraction` of that reference. Excluding the own column means the
vertical streak above/below a peak never counts against it. Only genuine
left/right neighbors do. It works in integer `(scan, TOF index)` space (no
calibration) and keeps/drops native points (no smoothing). Disable with
`--no-halo`.

## Experimental neighbor support

`--ms1-neighbor-radius`, `--prm-neighbor-radius`, and `--dia-neighbor-radius`
control temporal evidence independently. A radius of 1 uses up to one preceding
and one following **compatible observation**, plus the current observation.
`--neighbor-max-rt-gap` limits each neighbor's distance to the current frame
(default 5 seconds). No output points or intensities are imported from neighbors.

The legacy `--frame-half-width` CLI flag and `frame_half_width` TOML key remain
accepted for MS1. New recipes export `ms1_neighbor_radius`; conflicting TOML values
are rejected. Earlier running-average documentation described obsolete behavior:
the current pipeline sums for a keep mask and preserves native output intensities.

PRM/DIA event boundaries, geometry, acquisition settings, and calibration limits
are enforced. Larger radii change intensity-threshold and halo behavior; accuracy
and peak-boundary effects need downstream validation. See [full semantics](neighbors.md).

## Config file

Instead of (or alongside) flags, parameters can come from a TOML file:

```bash
dnoise <INPUT.d> <OUTPUT.d> --config dnoise.toml
```

Every key is optional and uses the same name as the flag with underscores. See
[dnoise.toml](../dnoise.toml) for a fully-commented example. Precedence is
**explicit CLI flag > config file > built-in default**, so a config sets the
baseline and individual flags override it for one run. Unknown keys are
rejected to catch typos.

```toml
mz_half_width = 3
min_feature_length = 7
iterations = 2
all_frames = false
# threads = 8
```

## Logging

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
adds trace, and `-q` limits output to warnings and errors. For fine-grained control set
`RUST_LOG` (e.g. `RUST_LOG=dnoise=debug`), which overrides the flags. When stderr is
an interactive terminal a progress bar is shown instead of periodic progress lines.

## Validating output

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

## Compression types

dnoise reads **compression type 2** input and always writes
type 2.

## Processing records and batch commands

Every successful output carries its own [processing history and reusable recipe](provenance.md).
Use `dnoise metadata output.d` to inspect history and `dnoise validate output.d`
to check SQLite metadata and every type-2 frame. An independent SDK check is
available with `python examples/validate_sdk.py --sdk /path/to/libtimsdata.so output.d`.
For unattended processing, see [batch manifests](batch.md).

Full validation uses up to four decoder tasks in the current Rayon pool, with
bounded batches read in physical file order. Normal runs use the pool selected
by `--threads`; `dnoise validate` uses Rayon's default pool (configurable with
`RAYON_NUM_THREADS`). Validation concurrency is capped at four to limit memory.
See [validation performance](validation-performance.md) for measured tradeoffs.

### Optional full validation

Full input/output decoding checks are **enabled by default**. To skip these
additional passes, use:

```sh
dnoise input.d output.d --skip-validation
dnoise batch jobs.json --skip-validation
```

Set `skip_validation = true` in TOML or a batch job's `config` for the same
behavior. `--validate` restores full validation and overrides configuration;
the two CLI flags cannot be combined. Batch flags override each job's setting.
The GUI has a **Full input/output validation** checkbox in advanced settings.

Skipping full validation still checks SQLite integrity, frame headers, offsets,
overlaps, parameters, and input/output paths. Processing still decodes the frames
it needs and uses checked encoding. Outputs are staged and structurally checked
before installation. Separate payload/count verification is skipped on both
input and output, including with `--in-place`; some corruption can go undetected.
Dry runs skip only the additional input decoding pass because they write no output.

Reports and processing history record `validation` as `full` or `structural`, and
saved recipes retain the setting. Older history without this field has no explicit
validation-mode record. `dnoise validate folder.d` always runs the full check,
regardless of any recipe stored in that folder.

## Synchro-, midia- and Slice-PASEF

These type-2 DIA-family acquisitions use the existing DIA options, including
`--denoise-msms` and `--dia-neighbor-radius`. Continuous one-scan window steps are
handled as scanning regions; static Slice/DIA windows stay separate. Neighbor
matching compares full trajectories. See [examples and limitations](acquisition-examples.md).
