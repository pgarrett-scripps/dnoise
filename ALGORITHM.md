# dnoise algorithm

This document describes the current Rust implementation. The source under
[`src/`](src/) is authoritative; [`docs/reference.md`](docs/reference.md) lists
the complete CLI and configuration interface.

## Data model

dnoise works directly in the integer coordinates stored by Bruker timsTOF
files. A decoded frame is represented by parallel arrays of ion-mobility scan,
TOF index, and intensity. Filtering in this native coordinate system avoids
floating-point binning. Calibration is used only by options expressed in
physical units, such as acquisition gates and region-of-interest crops.

The output directory retains the input metadata and tables, but
`analysis.tdf_bin` is rewritten with Bruker compression type 2. Frame offsets
and peak-count/intensity summaries in `analysis.tdf` are updated to match the
new binary.

## Default pipeline

The CLI defaults used in the benchmark are:

1. Filter MS1 frames with the iterative ion-mobility streak filter.
2. Apply the m/z-halo filter to the surviving MS1 points.
3. Apply the acquisition-appropriate MS1 gate when its geometry is present:
   the PASEF selection polygon for ddaPASEF, or the union of padded isolation
   windows for diaPASEF.
4. Copy MS/MS frames unchanged.
5. Encode every frame in its original order and update the database metadata.

MS/MS filtering, frame averaging, smoothing, centroiding, and cropping are
optional. The library API uses [`Stages::default`](src/params.rs), which leaves
all optional stages off; the CLI deliberately enables the benchmarked halo and
acquisition-aware MS1 gates.

## Ion-mobility streak filter

The core implementation is [`filter_once`](src/filter.rs), repeated by
[`filter_iterated`](src/filter.rs).

For every unique TOF index `c`, one pass:

1. Selects points whose TOF index is within
   `[c - mz_half_width, c + mz_half_width]`.
2. Sums their intensity independently in each ion-mobility scan.
3. Marks a scan occupied when its positive sum meets
   `min_window_intensity`.
4. Joins neighboring occupied scans across at most `max_internal_gap` empty
   scans.
5. Keeps a joined run when it contains at least `min_feature_length` occupied
   scans and its summed intensity meets `min_feature_intensity`.
6. Keeps points at TOF index `c` whose scan lies inside a retained run.

Bridged empty scans join runs but do not count toward the minimum feature
length. Iteration applies the same pass to the previous pass's survivors and
composes the keep mask back into the original point order.

The benchmarked MS1 defaults are:

| Parameter | Value |
|---|---:|
| `mz_half_width` | 3 TOF indices |
| `min_feature_length` | 5 occupied scans |
| `max_internal_gap` | 2 empty scans |
| `min_window_intensity` | 0 |
| `min_feature_intensity` | 0 |
| `num_iterations` | 2 |

## m/z-halo filter

After the streak filter, the default CLI compares each surviving point with
the maximum intensity in a surrounding scan/TOF box, excluding the point's own
TOF column. A point below `peak_fraction` of that off-column maximum is removed.
Excluding the own column prevents a point's vertical ion-mobility streak from
counting against it.

The benchmarked defaults are a peak fraction of 0.15, a TOF-index half-width
of 80, and a scan half-width of 2. `--no-halo` disables this stage.

## Acquisition-aware gates

The gates use acquisition geometry already stored in `analysis.tdf`:

- **ddaPASEF/PASEF MS1:** the padded IMS selection polygon removes survey
  points from precursor space the method never selects.
- **diaPASEF MS1:** the padded union of DIA isolation windows removes survey
  points from precursor space the method never fragments.
- **diaPASEF MS/MS:** when fragment frames are filtered, an isolation-window
  gate removes out-of-window scans and per-window filtering is enabled by
  default so unrelated windows are not linked.
- **ddaPASEF MS/MS:** an isolation-event gate enforces the recorded PASEF scan
  intervals when fragment frames are filtered.

Each gate is a no-op when its defining geometry is absent. The default physical
padding for both MS1 gate types is 5 Da in m/z and 0.05 1/K0 in mobility.

## Optional stages

- **MS/MS denoising:** `--denoise-msms` uses separately tuned streak-filter
  parameters. ddaPASEF precursor scans are combined before filtering;
  diaPASEF windows are filtered independently by default. This changes searched
  spectra and must be validated by re-searching.
- **Frame averaging:** `--frame-half-width` averages aligned MS1 coordinates
  across neighboring MS1 frames before filtering. It is experimental and off by
  default.
- **Smoothing:** `--smooth` replaces survivor intensities with local box
  averages before centroiding. Coordinates do not move.
- **Watershed centroiding:** `--watershed` grows intensity-ordered groups and
  emits one intensity-weighted integer centroid per group.
- **Box centroiding:** `--box-centroid` greedily tiles survivors into small,
  non-transitive boxes and emits one intensity-weighted centroid per box. It is
  mutually exclusive with watershed centroiding.
- **Crop:** m/z, mobility, retention-time, and intensity limits may be combined
  with denoising or applied alone with `--crop-only`. Frames outside an RT crop
  are emitted empty rather than deleted, preserving frame identifiers and table
  relationships.

Optional smoothing and centroiding operate only on frames processed by the
filter. The processing order is smoothing, then one selected centroider.

## Determinism and complexity

Filter and gate decisions are deterministic integer operations. Frames are
processed in parallel but written in their original order. Within one frame,
the streak filter sorts points by TOF and performs window lookups with binary
search; runtime therefore depends mainly on the number of points, unique TOF
indices, and enabled iterations. Memory use is bounded by decoded frames,
per-frame masks and profiles, and the ordered writer queue rather than the full
run.

## Limitations

- Input decoding currently supports Bruker compression type 2, not type 3.
- Acquisition gates assume the geometry stored in the run metadata is correct.
- Gates expressed in physical units use the run-level calibration; files with
  multiple calibration segments receive a warning because boundaries may be
  slightly offset in other segments.
- The published benchmark covers one instrument and two gradient lengths.
  Other instruments and sample types should be evaluated with `--dry-run` and
  downstream validation before routine use.

Run `dnoise --help` or read [`docs/reference.md`](docs/reference.md) for every
parameter, precedence rule, and operational safeguard.
