# Synchro-, midia- and Slice-PASEF support

Three public `.d` examples are stored locally under `test-data/acquisitions/`.
Raw datasets and generated outputs are excluded from Git; the checked source
manifest is `test-data/acquisition-examples.json`.

| Method | Public source | Example |
| --- | --- | --- |
| synchro-PASEF | [PXD037766](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD037766) | Three-minute four-protein run, timsTOF Pro 2, software 3.0.21. |
| midia-PASEF | [PXD073779](https://proteomecentral.proteomexchange.org/?pxid=PXD073779) | Five-minute NSCtrl run, timsTOF Ultra 2, software 6.0.8. |
| Slice-PASEF | [PXD073779](https://proteomecentral.proteomexchange.org/?pxid=PXD073779) | Five-minute 1-frame Slice run, timsTOF Ultra 2, software 6.0.8. |

The synchro example is a 43,632,133-byte ZIP. The midia and Slice examples are the
first complete `.d` directories in `midia.rar` and `Slice.tar.gz`; only the needed
archive prefixes were transferred. The full archives are approximately 188 and
619 GiB and were not downloaded. The manifest records source URLs, original member
names, actual native file sizes and SHA-256 checksums. Archive member paths and
file types were checked before extraction.

The source [synchro publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC9868879/)
and [comparison study](https://www.nature.com/articles/s41467-026-75845-5) identify
the acquisition methods. The examples are not sufficient to distinguish every
branded PASEF variant from metadata alone.

## Usage

All three use `Frames.MsMsType = 9`, compression type 2, and the standard
`DiaFrameMsMsInfo`/`DiaFrameMsMsWindows` tables. They are supported through the
existing DIA path; no new acquisition-name flag is required. CLI classification
continues to report the common **diaPASEF** family.

```bash
# Default: denoise MS1, preserve all native MS/MS points.
dnoise example.d ms1-output.d

# Experimental MS/MS filtering with one matching observation before and after.
dnoise example.d filtered-output.d --denoise-msms --dia-neighbor-radius 1
```

The existing `--msms-*`, halo, MS1 neighbor, crop, and optional postprocessing
settings apply. Fragment denoising and temporal support remain opt-in. Full
input/output validation remains enabled by default.

## Geometry handling

The inspected files encode their geometry differently:

- **Synchro:** four groups, 927 one-scan windows per group, with smoothly changing
  isolation m/z and collision energy.
- **midia:** ten groups, 944 one-scan windows per group. Groups encode different
  overlapping trajectories and must remain separate during neighbor matching.
- **Slice:** one group with 15 touching windows spanning 32–34 scans each. Their
  m/z centers and widths differ; touching scan endpoints do not make one window.

The checked DIA reader keeps static windows separate and combines successive
one-scan steps into scanning regions only when they are contiguous, move
monotonically in m/z, have equal widths (relative numerical tolerance `1e-9`),
and their isolation bands overlap. Gaps, direction reversals, width changes, and
nonoverlapping jumps end a region. This is a geometry-based interpretation of
recorded trajectories, not reconstruction of precursor-fragment identities.

Filtering and halo use each region's complete mobility extent. Treating every
one-scan metadata row as an independent event would otherwise reject ordinary
multi-scan features, including when using neighbors. Static Slice/DIA windows
remain independent even when they touch. This also corrects the prior file-reader
behavior that merged touching static intervals before per-window filtering.

Temporal support requires an exact match of the **entire recorded trajectory**,
including every scan interval, m/z center, width and collision energy. Matching
only the first window or bounding rectangle is insufficient. Calibration,
acquisition-setting and RT limits from [neighbor support](neighbors.md) still
apply. Neighbor matching indexes window groups once instead of expanding every
one-scan row for every frame.

When per-window processing is active, optional smoothing/centroiding also stays
inside these regions. `--no-dia-per-window` retains the explicit whole-frame
fallback when temporal support is off; a positive DIA neighbor radius forces
region boundaries. The original trajectory tables and frame times are copied
unchanged. Reports expose `stats.active_gates.dia_scan_varying` when scanning
regions are used for fragment filtering.

DIA metadata is validated before opening the frame reader, including interval
bounds, finite physical geometry, nonoverlap, group/frame references, and coverage
of nonempty DIA frames. This prevents malformed frame references reaching
unchecked indexing in the underlying reader.

See [repeatable validation](acquisition-validation.md) for the automated processing
and SDK comparison command.

## Validation limits

Synthetic coverage tests continuous one-scan trajectories, static touching
windows, trajectory discontinuities, exact internal m/z/energy matching, native
point retention, streaming/writer/dry-run parity, metadata preservation, and
malformed references. Existing PRM and standard DIA tests remain applicable.

Real-file checks compare default MS1-only and experimental radius-1 outputs with
the input using the independent Bruker SDK, including every frame, native-point
and intensity retention, unchanged default fragment spectra, metadata tables,
frame identities, calibration references and timestamps.

All three examples passed those checks. The measured radius-1 results, using
inherited filter settings and no smoothing/centroiding, were:

| Method | Frames per file | Input MS/MS points | Retained MS/MS points | Retained summed MS/MS intensity |
| --- | ---: | ---: | ---: | ---: |
| synchro-PASEF | 1,704 | 2,656,690 | 875,917 | 41.54% |
| midia-PASEF | 4,452 | 23,912,990 | 7,746,446 | 41.42% |
| Slice-PASEF | 4,455 | 392,833,100 | 195,301,363 | 51.44% |

Default MS1-only outputs preserved every fragment point and its intensity.
Detailed counts, commands, SDK checksum and run statistics are recorded in
[`test-data/acquisition-verification.json`](../test-data/acquisition-verification.json).
The full workspace passed 166 tests (two ignored), Clippy with warnings denied,
and formatting checks; the root library/CLI also passed the Rust 1.85 check.

This establishes file compatibility and native-subset behavior for these three
examples. It does **not** validate peptide identifications, fragment ratios,
chromatographic shapes, quantitative accuracy, or deconvolution. The midia example
is an NSCtrl run, not a quantitative performance benchmark. Substantial summed
fragment intensity can be removed with the inherited filter settings. Do not
interpret file reduction or successful decoding as quantitative acceptance.
