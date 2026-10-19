# dnoise roadmap

Post-submission plan · October 2026 · Released version **0.6.0**.

**Steps 1 to 3 and both extensions from step 5 are implemented and released.**
Three releases have shipped since 0.3.0, each through the gate in
[RELEASE.md](RELEASE.md) and archived on Zenodo: 0.4.0 (2026-09-23), 0.5.0
(2026-09-25) and 0.6.0 (2026-10-11). The manuscript is pinned to the v0.4.0 tag
and its benchmark settings are unchanged, so no later release moves a committed
manuscript number. What is still open is the manual GUI smoke checks in step 3
and the low-input research track in step 4.

Defaults changed twice since 0.3.0 and the workspace gained a crate, so read
[CHANGELOG.md](CHANGELOG.md) before comparing output across versions.

| Release | What changed for users |
|---|---|
| 0.4.0 | MS1 acquisition gates decide per feature instead of per point, and their fixed pads are gone |
| 0.5.0 | New `dnoise-core` crate, halo `peak_fraction` 0.10, MS1 gate pads 3.0 Th and 0.015 1/K0, MS1 stage order streak, gate, halo |
| 0.6.0 | Frames read with `sage-plus-tdf` instead of `timsrust-tdf`, minimum Rust 1.97 for the library and CLI |

## 1. Protect files · implemented

The original `--force` path could delete input when input and output matched.
The CLI, GUI, and high-level library writer now share a protected output lifecycle.

- Reject equal, overlapping, and aliased input/output paths.
- Process and validate in an owned temporary folder before replacing output.
- Preserve existing files on errors and cancellation, and restore a backup after
  a failed installation or report its recovery location.
- Exercise overwrite, cancellation, copy failure, and failed-swap regressions.

**Evidence:** the safety tests pass. Default output is deliberately no longer
byte-identical to the original submission, because 0.4.0 and 0.5.0 changed the
MS1 defaults, so reproducing the paper means building the v0.4.0 tag. Within a
version the output path is stable: the 0.5.0 core split and the 0.6.0 reader swap
were both checked to leave output unchanged. Replacement is recoverable but is
not an atomic directory exchange, so a crash can still require manual recovery.

## 2. Validation, provenance, and maintenance · implemented

Completed outputs now carry a processing record inside the `.d` folder.

- Generate small DDA/DIA fixtures for read/filter/write/reopen and streaming
  parity tests, and reject malformed frames through checked decoding.
- Provide `dnoise validate` and a separate Bruker SDK validation script. Full
  validation is on by default, and `--skip-validation` retains structural checks.
- Resolve CLI/GUI settings through shared code. Write `dnoise.provenance.json`
  and `dnoise.config.toml`, including history, parameters, versions, and statistics.
  Since 0.5.0 the record also carries `unit_equivalents`, the physical reading of
  each raw-unit parameter under that run's calibration.
- Update affected dependencies and the GUI framework. Configure explicit compiler
  selection, platform tests, dependency audits, and monthly update PRs.

**Evidence:** 210 workspace tests pass across 20 test binaries, with 2 ignored
because both need a real `.d` folder: the streaming parity test, which runs when
`DNOISE_PARITY_DOTD` points at one, and the `RunContext` doc example.
`dnoise-core` still checks on Rust 1.85, the library and CLI on 1.97, and the GUI
on 1.95, each enforced by its own CI job. The audit job runs `cargo audit` with no ignored advisories: 0.6.0's move
to `sage-plus-tdf` dropped the unused cloud-storage stack that RUSTSEC-2026-0194
and RUSTSEC-2026-0195 came in through. No custom TDF tables or metadata keys are
added.

## 3. Easier daily use · implemented

The desktop workflow supports selection, estimation, processing, and inspection.

- Native folder pickers and visible input/configuration errors and warnings.
- Batch collision checks, cancellation, and retry of unfinished files.
- Completed raw/output frame comparisons using matching axes and color scales.
- Configurable frame batching, with the original 2,048-frame default retained.
  Smaller batches reduced memory in a local probe but took longer.

CI now runs formatting, clippy, and the workspace tests on Linux, macOS, and
Windows for every push and pull request, which covers the GUI build on all three.
Local headless UI tests pass. Full input and output validation uses Rayon: local
complete runs took 23 to 31% less time than serial validation, with identical
output. See the [measurements](docs/validation-performance.md). Cancellation and
performance timings are not guarantees for every acquisition or disk.

**Still open:** native GUI and dialog smoke checks on macOS and Windows are
manual, so nothing in CI catches a picker or dialog regression on those two
platforms. A Windows installer with a `.d` folder right-click entry exists on the
unmerged `worktree-windows-installer` branch, one commit ahead of main and 40
behind it, and needs a test on real Windows 11 before it lands. Runtime and memory measurements
still cover local probes rather than a representative set of runs.

## 4. Better low-input performance · deferred research

The manuscript identifies weak-signal losses in sparse and low-input runs.
Use those findings to guide a separately validated experiment.

- Compare relaxed thresholds and individual stages against raw data and paper defaults.
- Address calibration-segment handling for physical-unit gates and crops.
- Evaluate held-out, replicated runs using identification, quantification,
  missingness, and retention by abundance alongside file-size reduction.
- Offer an opt-in low-input preset only after it meets predefined fidelity goals.

A relaxed streak setting was evaluated on the manuscript's own low-input runs
during the revision, and no preset was added to the code. Those runs were used to
describe the behavior, not held out from it, so the decision gate below still
stands: demonstrate a useful fidelity/reduction tradeoff on data that were not
used to tune the settings. Allow 4 to 8 calendar weeks for an initial study,
depending on data and compute. This work is not part of the current implementation.

## 5. Integrations · batch workflows and embedding implemented

The first extension is a portable JSON batch manifest for unattended lab runs.
The GUI can export a queue, and `dnoise batch` preflights writing runs, executes
jobs, records per-file results, and exits nonzero if any job fails.
[Workflow examples](docs/batch.md) cover relative paths, reports, and retries.

The second landed in 0.5.0. `dnoise-core` carries the in-memory stages with no
I/O and no native dependencies, so another program can read frames its own way
and call `dnoise_core::Denoiser` per frame. CI fails if the crate gains an I/O or
native dependency. sage-plus is the motivating consumer. In-process use that
still wants dnoise to write the output goes through `RunContext` in `dnoise`.

Python bindings remain contingent on demand beyond those two paths. Additional
formats require real unsupported inputs and independent validation tools.
Algorithm or GPU work needs a measured bottleneck or scientific justification.

## Ongoing maintenance

| When | Action |
|---|---|
| Weekly | Triage bugs and reviewer requests, about 20 minutes |
| Monthly | Review dependency updates and advisories, about 60 to 90 minutes |
| Every release | Run platform/SDK checks, record hashes, update changelog, archive |
| On acceptance | Add the journal DOI and final citation |

[Implementation and validation evidence](docs/implementation-validation-2026-09-09.md)
· [Original code review](docs/maintenance-review-2026-09-09.md)
· [Metadata and recovery](docs/provenance.md).
