# Targeted proteomics: prm-PASEF

The current development version supports **compression type-2 prm-PASEF `.d`
folders**, with MS1-only processing by default and opt-in experimental fragment
denoising:

```bash
dnoise input.d output.d
dnoise input.d output.d --dry-run
dnoise input.d output.d --denoise-msms
dnoise validate output.d
```

The CLI, GUI, streaming library, and provenance reports recognize PRM. MS1 uses
the existing filter and optional stages. PRM fragment coordinates and intensities
pass through unchanged unless MS/MS denoising or an explicit crop is requested. Fragment records
are re-encoded, so their compressed bytes need not match. Original target,
isolation-event, scheduling, calibration, and auxiliary tables are retained;
only the usual frame offsets and peak/intensity summaries are updated.

## Acquisition policy

| Acquisition | MS1 filtering | Discovery gates | MS/MS denoising |
|---|---|---|---|
| ddaPASEF / diaPASEF | Existing behavior | Existing behavior | Existing optional paths |
| prm-PASEF | Supported | Disabled | Experimental, opt-in per event |
| Mixed / unknown | Supported at the frame level | Disabled | Rejected |

Mixed means more than one nonzero `Frames.MsMsType`, including runs that contain
PRM alongside DDA or DIA. A conservative policy applies to the whole mixed run.
The file writer and streaming context share this policy. PRM `--denoise-msms`
uses the recorded isolation events; `--all-frames` remains rejected to avoid
applying MS1 parameters to fragments. Mixed/unknown runs reject both options.
Crop-only processing does not run denoising. Ordinary crop
bounds also apply to MS/MS and can change targeted quantitative results.

Files without MS1 frames report that there is no MS1 denoising opportunity.
Do not estimate savings from frame counts: measure bytes or run the tool.

Detection checks PRM event frame/target references, finite positive isolation
m/z and width, nonnegative collision energy, and bounded, nonoverlapping scan
intervals with exclusive ends. Adjacent events remain distinct. Nonempty PRM
frames require events; empty frames may lack them. Conflicting DDA/DIA events on
PRM frames and invalid optional measurement-mode frame references are rejected.
Empty PRM tables in discovery files do not activate PRM handling. These metadata
checks also run during `validate`, streaming setup, and `--skip-validation`.

## Experimental PRM MS/MS filtering

`--denoise-msms` runs the streak filter and enabled halo independently in each
recorded `[ScanNumBegin, ScanNumEnd)` interval. Adjacent events are never merged,
and repeated observations stay independent by default. Opt-in
`--prm-neighbor-radius N` supplies evidence from nearby matching target observations
for the keep decision only; see [neighbor support](neighbors.md).
Points outside all recorded intervals are dropped. MS1 processing is unchanged.

The existing `--msms-*` settings apply: defaults are a TOF half-width of 3,
minimum 3 occupied scans, maximum internal gap of 8, one iteration, and zero
intensity floors. These are inherited settings, not a PRM-optimized preset.
For example, `--msms-min-feature-length 2` changes the minimum for each event.

Surviving points keep their native coordinates, intensities, and frame times.
If optional smoothing or centroiding is enabled, it also runs independently
inside each event; those stages intentionally change intensities or coordinates.
Discovery gates and `dia_per_window` settings cannot override PRM event boundaries.
Crop-only runs bypass all denoising, including PRM event filtering.

CLI/GUI reports include an experimental warning, and
`stats.active_gates.prm_per_event` records whether the path was enabled. Validate
extracted fragment areas, ratios, chromatographic shape, missingness, and replicate
precision before routine use. Structural compatibility and retention of native
survivor intensities do not establish quantitative fidelity.

## Default MS1-only validation: September 9, 2026

One unmodified public file was downloaded from
[PXD049405](https://www.ebi.ac.uk/pride/archive/projects/PXD049405):
`200ngHeLa_top2000_apop_PRM_45m_6_Slot1-2_11-1-2023_7479.d.zip` (458,282,811 bytes).
It contains 82 targets and 64,048 isolation events, with `MsMsType=10` for PRM.
The original SQLite schema supplied the synthetic fixture's PRM table layout.

The release-profile build ran the ordinary CLI defaults with four threads and
full input/output validation. A separate Python/ctypes comparison used the
locally available Bruker SDK to decode **every frame of both files**, check scan
and peak counts, and compare every PRM scan's TOF/intensity arrays exactly.

| Measurement | Result |
|---|---:|
| Frames per input/output | 33,088 |
| MS1 frames | 87 |
| PRM frames with identical SDK-decoded data | 33,001 |
| PRM points preserved exactly | 112,585,040 |
| Input binary | 437,542,912 bytes |
| Output binary | 400,884,828 bytes |
| Binary reduction | 8.38% |
| Processing including full validation, local four-thread run | 8.67 seconds |

All 18 non-`Frames` tables and the SQLite schema matched between input/output.
Frame IDs, retention times, acquisition modes, scan counts, accumulation times,
and calibration references matched. Input `analysis.tdf` and `analysis.tdf_bin`
still matched the original archive CRCs after processing.

Reading actual frame headers, MS1 records occupied 47,083,145 bytes (10.76% of
the input binary); PRM records occupied 388,321,781 bytes. The remainder was
header/padding space. Overall reduction includes filtering, re-encoding, and
layout changes; it is not exclusively deleted MS1 signal.

SDK library SHA-256 used for this comparison:
`a97086132c64bed75496f7da02aa43a741a8dba4b1dc031e194cd9c51a6f3e9d`.
For independent frame-count/decoding checks, the existing helper accepts a local
SDK library:

```bash
python examples/validate_sdk.py --sdk /path/to/libtimsdata.so input.d output.d
```

That helper checks readability and counts; the exact PRM-array and SQL-table
comparisons above were additional checks in this validation run.

## Experimental MS/MS validation: September 9, 2026

The same public run was processed with `--denoise-msms --threads 4`, retaining
the inherited MS/MS defaults and the default halo. Full input/output validation
was enabled; no smoothing, centroiding, or cropping was requested.

| Measurement | Result |
|---|---:|
| Output binary | 121,316,380 bytes |
| Binary reduction versus raw | 72.27% |
| Retained PRM points | 31,625,315 of 112,585,040 |
| Total PRM intensity retained | 41.28% |
| Processing including full validation | 21.93 seconds |

The Bruker SDK independently read all 33,088 frames from both files. For every
PRM frame, the output was checked as a subset of the input's `(scan, TOF,
intensity)` points, and every retained point lay inside a recorded event. All
18 non-`Frames` tables, the SQLite schema, and frame identity/time/calibration
metadata were unchanged. The input still matched the archive CRCs.

**The initial settings removed 58.72% of summed PRM intensity.** This sum includes
background and all measured fragments; it is not a measure of peptide-level
quantitative bias. Its substantial change means these settings should not be
treated as a validated PRM preset. Peptide/transition-level peak areas and
replicate precision remain untested. The large file reduction alone is not an
acceptance criterion for quantitative use.

Synthetic tests cover default fragment and metadata preservation, writer/streaming
and dry-run parity, touching targets, sparse target IDs, empty frames/tables,
mixed runs, explicit crops, malformed metadata, all-frames rejection, CLI behavior,
and preservation of an existing output on failure. PRM MS/MS tests also verify
that targets are not pooled across frames by default and halo/smoothing/centroiding do not
mix adjacent events. Workspace tests (152 passed), Clippy,
library-only compilation, and the Rust 1.85 all-target check passed.

This establishes compatibility for one acquisition/software example. Skyline or
SpectroDive extraction, other instruments/software versions, and quantitative
effects of filtering have not been evaluated. Conventional non-mobility PRM,
SRM/MRM, and compression type 3 remain unsupported.


Neighbor support is now available via `--prm-neighbor-radius`; the radius-1
comparison and current checks are documented in [neighbor validation](neighbors.md#validation-run-2026-09-09).
The original radius-0 measurements above remain the baseline.
