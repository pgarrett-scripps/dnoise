# dnoise roadmap

Post-submission plan · September 2026 · Development version **0.3.0**.

**Steps 1–3 and the batch-workflow extension from step 5 are implemented locally.**
Release verification remains open. The paper’s v0.1.0 tag and benchmark settings
are unchanged; the low-input research track is deferred.

## 1. Protect files — implemented

The original `--force` path could delete input when input and output matched.
The CLI, GUI, and high-level library writer now share a protected output lifecycle.

- Reject equal, overlapping, and aliased input/output paths.
- Process and validate in an owned temporary folder before replacing output.
- Preserve existing files on errors and cancellation; restore a backup after a
  failed installation or report its recovery location.
- Exercise overwrite, cancellation, copy failure, and failed-swap regressions.

**Evidence:** safety tests pass, and default DDA/DIA output binaries match the
paper release byte for byte. Replacement is recoverable but is not an atomic
directory exchange; crashes can require manual recovery.

## 2. Validation, provenance, and maintenance — implemented

Completed outputs now carry a processing record inside the `.d` folder.

- Generate small DDA/DIA fixtures for read/filter/write/reopen and streaming
  parity tests; reject malformed frames through checked decoding.
- Provide `dnoise validate` and a separate Bruker SDK validation script. Full
  validation is on by default; `--skip-validation` retains structural checks.
- Resolve CLI/GUI settings through shared code. Write `dnoise.provenance.json`
  and `dnoise.config.toml`, including history, parameters, versions, and statistics.
- Update affected dependencies and the GUI framework. Configure explicit compiler
  selection, platform tests, dependency audits, and monthly update PRs.

**Evidence:** 142 workspace tests pass; the core still checks on Rust 1.85.
The current dependency audit has zero known vulnerabilities and one transitive
unmaintained-package warning. No custom TDF tables or metadata keys are added.

## 3. Easier daily use — implemented; release checks open

The desktop workflow now supports selection, estimation, processing, and inspection.

- Native folder pickers and visible input/configuration errors and warnings.
- Batch collision checks, cancellation, and retry of unfinished files.
- Completed raw/output frame comparisons using matching axes and color scales.
- Configurable frame batching, with the original 2,048-frame default retained.
  Smaller batches reduced memory in a local probe but took longer.

**Before release:** run native GUI/dialog smoke checks on Linux, macOS, and Windows;
exercise the configured platform CI; expand runtime/memory and cancellation
measurements to representative runs. Local headless UI tests pass. Full input
and output validation now uses Rayon: local complete runs took 23–31% less time
than serial validation, with identical output. See the [measurements](docs/validation-performance.md).
Cancellation and performance timings are not guarantees for every acquisition or disk.
Improve installers if launch problems recur.

## 4. Better low-input performance — deferred research

The manuscript identifies weak-signal losses in sparse and low-input runs.
Use those findings to guide a separately validated experiment.

- Compare relaxed thresholds and individual stages against raw data and paper defaults.
- Address calibration-segment handling for physical-unit gates and crops.
- Evaluate held-out, replicated runs using identification, quantification,
  missingness, and retention by abundance alongside file-size reduction.
- Offer an opt-in low-input preset only after it meets predefined fidelity goals.

**Decision gate:** demonstrate a useful fidelity/reduction tradeoff on data that
were not used to tune the settings. Allow 4–8 calendar weeks for an initial study,
depending on data and compute. This work is not part of the current implementation.

## 5. Integrations — batch workflows implemented

The first extension is a portable JSON batch manifest for unattended lab runs.
The GUI can export a queue; `dnoise batch` preflights writing runs, executes jobs,
records per-file results, and exits nonzero if any job fails.

[Workflow examples](docs/batch.md) cover relative paths, reports, and retries.
Python bindings remain contingent on an in-process analysis workflow; additional
formats require real unsupported inputs and independent validation tools.
Algorithm or GPU work needs a measured bottleneck or scientific justification.

## Ongoing maintenance

| When | Action |
|---|---|
| Weekly | Triage bugs and reviewer requests; about 20 minutes |
| Monthly | Review dependency updates and advisories; about 60–90 minutes |
| Every release | Run platform/SDK checks, record hashes, update changelog, archive |
| On acceptance | Add the journal DOI and final citation |

[Implementation and validation evidence](docs/implementation-validation-2026-09-09.md)
· [Original code review](docs/maintenance-review-2026-09-09.md)
· [Metadata and recovery](docs/provenance.md).
