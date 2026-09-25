# Processing metadata

Every successful output contains two files alongside `analysis.tdf` and
`analysis.tdf_bin`:

| File | Purpose |
|---|---|
| `dnoise.provenance.json` | Completion marker and processing history |
| `dnoise.config.toml` | Exact executable recipe for the latest operation |

The JSON records schema/software versions, build revision, completion time,
input basename, acquisition type, all active stage parameters, active gates,
crop settings, point/intensity statistics, actual binary sizes, and warnings.
It labels crop-only operations separately from denoising.
`unit_equivalents` lists each active raw-unit parameter (TOF indices, scans)
with what it meant on this run's calibration: TOF half-widths in ppm at m/z 400,
800 and 1200, scan counts in 1/K0 at the mid scan. The run log prints the same
values, one `units:` line per parameter.
New records also identify the validation mode (`full` or `structural`). Structural
mode means `--skip-validation` was used; completion is not a claim that every
input/output payload passed a separate decoding check. The TOML retains this choice.

The TOML contains resolved values, including the effective integer window when
ppm conversion was requested. It reproduces the recipe on the same input; moving
that integer window to another instrument does not preserve a fixed ppm tolerance.
Library-only builds without the `config` feature still write the complete JSON,
but do not export TOML or retain a stale TOML recipe copied from the input.

```sh
dnoise metadata processed.d
dnoise original.d reproduced.d --config processed.d/dnoise.config.toml
```

Reprocessing appends a history entry and replaces the latest TOML recipe.
CLI/GUI runs warn when input already contains history. Dry runs write neither
file, although an explicitly requested external report may be written.
Malformed or unsupported existing history is rejected rather than discarded.

**A config file alone is not a completion marker.** Both sidecars are written
inside a temporary output after processing and validation, then installed with
the completed acquisition. Failed or cancelled runs do not mark the input.

## Why a sidecar instead of custom TDF fields?

This release adds no custom tables or keys to Bruker's database. It changes only
the normal frame offsets/counts and compression metadata required by the writer.
The sidecars are readable without SQLite and can evolve independently of the
vendor schema. An output with the sidecars was read successfully by the Bruker
SDK during development; this is not a guarantee for every downstream tool.

Keep both sidecars when copying or depositing the folder. Tools that copy only
the two vendor files will lose the processing record. Absence of a marker does
not establish that data are raw: older dnoise versions wrote no marker, and
sidecars can be removed. The marker is provenance, not tamper-proof certification.

## File safety and recovery

Normal input/output paths must be disjoint, including symlink aliases and
parent/child directories. Use explicit `--in-place` to replace an input.
The output's parent directory must already exist. Symlinks inside acquisitions
are rejected rather than followed while copying.

Each writer uses an exclusive `<output>.dnoise-lock` and an owned temporary
sibling. Existing output is moved to a unique backup only when installation is
ready. A failed install restores that backup; if restoration fails, the error
identifies where the original remains. A crash can leave staging/backup folders
and a lock: inspect them before manually recovering or removing anything.
This is recoverable replacement, not a single atomic directory exchange or a
guarantee against power loss or unrelated programs changing files concurrently.

External `--report` paths must not exist or overlap an acquisition. Reports now
have `schema_version = 1`; consumers of older unversioned reports should check
the new schema. Output provenance remains available if an external report fails.

## Intensity and neighbor QC

Reports declare `intensity_basis: "stored_integer"`; these are the integers in
the frame payload, before any normalization a vendor reader may apply.
Reports now include `raw_ms1_summed_intensity`, `kept_ms1_summed_intensity`,
`raw_msms_summed_intensity`, and `kept_msms_summed_intensity` under `stats`.
`ms1_intensity_retained_pct` and `msms_intensity_retained_pct` divide output by
input intensity, returning JSON `null` when the input total is zero. Sums cover
all processed frames, including intensity removed by an RT crop; dry-run sampling
covers only sampled central frames. RT-excluded nonempty frames are now decoded
to obtain those totals, which can add work to cropped runs. Smoothing/centroiding,
if explicitly enabled, can alter output intensities; these statistics describe
actual output, not an assurance of native-intensity preservation.

`ms1_neighbor_usage` and `msms_neighbor_usage` each report:

| Field | Meaning |
| --- | --- |
| `events` | Nonempty central events evaluated with temporal support enabled; one event per MS1 frame |
| `events_without_neighbors` | Evaluated events that received no matching neighbor points |
| `neighbors_used` | Total other observations supplying points inside matching events after RT, empty-frame and crop exclusions |
| `max_neighbors_used` | Largest count of such observations for one central event |

The central observation is excluded from neighbor counts. Reusing one neighbor
for different central events counts each use; these are not unique frame counts.
Empty central frames/events, crop-only runs and radius-zero paths contribute no
events. Counts measure supplied evidence, not causal effects on retained peaks.

## Calibration restrictions

Physical DIA MS1 windows and DDA selection polygons require a single referenced
m/z and mobility calibration. If multiple references are present, an applicable
gate now errors instead of using a run-level approximation. Disable that gate
explicitly (`--no-dia-ms1-window` or `--no-ms1-polygon`) to use raw-coordinate
filtering. Gates with no applicable geometry remain inactive without an error.

Physical crops reject multiple calibrations only in the requested dimension:
m/z bounds require one m/z calibration, mobility bounds require one mobility
calibration. RT-only and intensity-only crops remain available. Crop-only mode
does not build acquisition gates. The same checks protect writer, dry-run and
streaming paths; `--skip-validation` does not bypass them. Temporal matching
continues to separate calibration changes.
