# Experimental neighbor support

Neighbor support sums nearby compatible observations to decide which native
points survive. It does **not** write a summed or averaged spectrum, add peaks
from another frame, or change retained intensities. Explicit smoothing or
centroiding remains a separate operation that can change output points.

```bash
dnoise input.d output.d --denoise-msms \
  --prm-neighbor-radius 1 --dia-neighbor-radius 1 --ms1-neighbor-radius 1 \
  --neighbor-max-rt-gap 5
```

| Option / TOML key | Default | Meaning |
| --- | --- | --- |
| `--ms1-neighbor-radius` / `ms1_neighbor_radius` | 0 | Previous/next compatible MS1 observations per side. |
| `--prm-neighbor-radius` / `prm_neighbor_radius` | 0 | Previous/next compatible observations of the same PRM target per side. |
| `--dia-neighbor-radius` / `dia_neighbor_radius` | 0 | Previous/next compatible DIA isolation-window observations per side. |
| `--neighbor-max-rt-gap` / `neighbor_max_rt_gap` | 5 | Maximum absolute retention-time distance from the current frame, in **seconds**; finite and positive. |

Radius 1 uses up to three observations: previous, current, next. Radius 2 uses up
to five. PRM targets and DIA windows may recur many physical frames apart; unrelated
frames do not consume slots. Empty matching observations do consume slots but
supply no points. At run edges, or when observations exceed the RT limit, the
window shrinks without backfilling from farther away.

PRM requires `--denoise-msms`; DIA requires `--denoise-msms` or `--all-frames`.
A positive radius with fragment filtering disabled is rejected on the matching
acquisition. PRM/DIA radii are ignored on other acquisition types, allowing a shared
batch recipe. DDA retains its existing precursor-based pooling.

## Matching and boundaries

- PRM requires the same target ID, absolute scan interval, isolation m/z and width,
  and collision energy. DIA requires the same exact interval and isolation
  geometry/energy; equivalent windows can match across window-group IDs. For
  synchro/midia scanning regions the entire scan-dependent trajectory must match;
  see [acquisition examples](acquisition-examples.md).
- Frame scan counts and acquisition settings must match: scan mode, polarity,
  accumulation time and ramp time when present, plus PRM measurement mode when
  recorded. No interpolation, mobility shifting or m/z recalibration is performed.
- Both calibration references must be present and non-null. Any change in the
  MzCalibration/TimsCalibration references starts a new segment; neighborhoods
  cannot cross that boundary, even if a later frame returns to earlier references.
- Frame retention times must be finite and nondecreasing. Each neighbor must be
  within the configured RT distance of the current frame. RT-cropped frames cannot
  supply evidence; point-level crops still run after filtering. Crop-only mode
  bypasses neighbor processing entirely.
- PRM and DIA events remain distinct even when their scan intervals touch.
  Continuous one-scan quadrupole steps in synchro/midia form scanning regions;
  static Slice/DIA windows stay separate. Temporal support forces region-local
  DIA filtering despite `--no-dia-per-window`.
  Halo and optional smoothing/centroiding also stay inside each current event.
  Out-of-event fragment points are discarded when temporal support is enabled.
  Invalid, overlapping, missing or orphaned DIA event metadata is rejected.

## Filtering and quantitative limits

Within an event, native intensities at identical `(scan, TOF)` bins are summed,
then the existing streak filter and optional halo produce a set of surviving
coordinates. That set is applied **only to points already present in the current
frame**. The existing MS1 or MS/MS filter parameters still apply. As in the DDA
combiner, summed bins saturate at `u32::MAX` before filtering.

Intensity thresholds apply to the **sum**, without dividing by neighborhood size.
A larger radius therefore changes their effective stringency; run edges have less
support. With zero intensity thresholds, the benefit comes from combined scan
occupancy and halo evidence, not merely a higher summed intensity. Strong neighbor
signals can also change halo decisions, so retention need not increase monotonically.

This remains experimental. Native-intensity preservation does not guarantee
fragment-area, ratio, peak-shape or peptide-level quantitative preservation.
Neighbor evidence may preserve noise near peak boundaries. Compare radii 0, 1 and
2 downstream before choosing a production setting; the 5-second limit is a guard,
not a validated quantitative preset.

## Interfaces and verification

The CLI alias `--frame-half-width` and legacy TOML key `frame_half_width` remain
accepted for MS1. Conflicting old/new TOML values are rejected; explicit CLI flags
override both. New saved recipes use `ms1_neighbor_radius`. Library callers keep
`Stages::frame_half_width` for MS1 and set `Stages::neighbors: NeighborParams` for
PRM/DIA and the shared time limit. The GUI exposes all four controls in Advanced
settings. Reports include resolved settings and active-neighbor flags.

Metadata is indexed once. Each frame worker decodes only selected supporting
frames, reusing decoded neighbors across that frame's events. There is no global
cache of all spectra. Cost grows with radius and active workers. Streaming, file
writing and dry runs use the same decision path; sampled estimates may read
unsampled source frames for support.

Synthetic tests cover both-sided support, radii 0/1/2, native-point/intensity
preservation, summed thresholds, empty frames, run edges, exact target/window
matching, calibration and RT boundaries, crop-only behavior, event-local halo and
postprocessing, malformed DIA metadata, CLI aliases/overrides, replayable recipes,
and parity across streaming, dry runs, thread counts and writer batch boundaries.

## Validation run (2026-09-09)

The workspace suite passed **161 tests**, including the new synthetic coverage for
all three workflows. Clippy, library-only tests (75), Rust 1.85 library/CLI
all-target compilation, and formatting checks passed. The GUI builds and its
existing rendering tests passed on the workspace toolchain.

The public PXD049405 PRM run described in [targeted support](targeted.md) was
processed with all radii set to 1, a 5-second RT limit, four workers, default
streak/halo settings, full input/output validation, and no smoothing or centroiding.
The DIA radius was inactive because this was PRM data.

| Measurement | Radius 0 | Radius 1 |
| --- | ---: | ---: |
| Retained PRM points | 31,625,315 | 41,242,356 |
| Retained MS1 points | 4,634,747 | 6,272,691 |
| Output binary bytes | 121,316,380 | 159,015,500 |
| Binary reduction | 72.27% | 63.66% |
| Run time including validation (seconds) | 21.93 | 41.09 |

These measurements show greater signal retention at increased processing cost
and output size. They do not establish quantitative accuracy. DIA and MS1 have
synthetic coverage; the real dataset is a PRM acquisition containing MS1 frames.

Independent Bruker SDK verification read all **33,088 frames**
in both input and output and confirmed that retained points in every MS1 and PRM
frame kept their native coordinates and intensities. All retained PRM points stayed
inside their recorded events; all 18 non-Frames metadata tables,
SQLite schema, frame timing and calibration references were unchanged. The input
still matched the downloaded archive CRCs.

Using SDK-decoded intensities, PRM summed-intensity retention increased from
**41.28% to 48.38%**. About
**51.62%** was still removed; these totals include
background and all fragments and do not measure peptide-level quantitative bias.
