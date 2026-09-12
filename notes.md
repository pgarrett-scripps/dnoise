
Adds prm-PASEF and the scanning diaPASEF variants, a streaming library API, and
in-file provenance.

Denoised output is byte-identical to 0.1.0 for ddaPASEF on all three arms
(MS1-only, MS1 + MS/MS, watershed centroider) and for diaPASEF MS1-only,
verified by running both builds over the benchmark acquisitions and comparing
the stored frame data and metadata.

**diaPASEF MS1 + MS/MS output changed.** The DIA window-grouping rework, which
keeps static touching windows separate and groups continuous one-scan
trajectories, filters fragments slightly harder: roughly 0.37% fewer MS/MS
points are retained, across about half the MS/MS frames. MS1 is unaffected.
Anyone reproducing published diaPASEF fragment-denoising numbers from 0.1.0
should pin 0.1.0 or re-derive them.

### Fixed
- Write empty (0-peak) frames the long way, as a header plus a compressed
  all-zero scan table, rather than the header-only 8-byte record. timsrust
  zstd-decodes whatever follows the header unconditionally and fails on an
  absent payload, so a single emptied frame made the whole `.d` unreadable to
  every downstream consumer, Sage included. MS/MS denoising empties frames
  routinely. The reader still accepts both shapes.

### Added and changed
- Reject physical gates/crops that require unsupported multi-calibration conversions;
  retain raw-coordinate filtering and RT/intensity-only crops.
- Report separate MS1/MS/MS intensity retention and actual neighbor usage, including
  intensity removed by RT cropping in the input totals.
- Add a reusable real-acquisition comparison command with independent SDK checks,
  pinned input checksums, exact executable hashes, recipes, and failure reports.

- Add prm-PASEF support: validate PRM target/event metadata, detect PRM and mixed
  acquisitions in CLI/GUI/reports, and preserve fragment points by default.
  `--denoise-msms` enables experimental filtering within each PRM isolation
  event, including event-local halo/smoothing/centroiding. Cross-frame pooling
  is disabled by default. Record the active path in provenance. Disable discovery gates on
  PRM/mixed/unknown data; mixed/unknown runs reject fragment denoising. Validate a
  public PRM run with the Bruker SDK; see [targeted support](docs/targeted.md).
- Add experimental PRM/DIA/MS1 neighbor radii in CLI, TOML, library stages, GUI,
  and provenance. Sum compatible nearby observations for keep decisions without
  importing peaks or intensities. Enforce event boundaries, calibration segments,
  acquisition settings, and a configurable RT-distance limit (default 5 seconds).
  Keep `--frame-half-width` / `frame_half_width` as MS1 aliases; correct obsolete
  averaging documentation. See [neighbor support](docs/neighbors.md).
- Support public synchro-PASEF, midia-PASEF, and Slice-PASEF examples through the
  DIA path. Keep static touching windows separate, group continuous one-scan
  trajectories for filtering, and match complete trajectories for temporal
  support. Validate DIA references before opening the decoder. Record scanning
  filtering with `ActiveGates::dia_scan_varying`; see [examples](docs/acquisition-examples.md).
- Protect input/output paths and install only validated temporary outputs, with
  ownership-aware cancellation cleanup and recoverable in-place replacement.
- Store processing history and exact configuration inside every completed `.d`.
  Add shared configuration resolution, versioned reports, checked codec APIs,
  synthetic DDA/DIA integration tests, and `validate`/`metadata` commands.
- Add native folder pickers, batch collision checks, retry controls, completed
  frame comparisons, and per-frame cancellation checks. Expose frame batch size.
- Add portable JSON batch manifests and `dnoise batch`, including per-file results.
- Upgrade affected dependencies and the GUI framework. Add dependency audits,
  monthly update PRs, explicit compiler selection, and platform test gates.
- Parallelize full-frame validation with Rayon using bounded read batches and
  at most four decoder tasks; retain structural, content, and cancellation checks.
- Make full input/output validation optional via `--skip-validation`, shared
  configuration, and a GUI checkbox. Keep full validation on by default, retain
  structural/file-safety checks, and record the mode in provenance and recipes.

Library migration: `Stages` adds `neighbors: NeighborParams`; `ActiveGates` adds
`ms1_neighbors`, `prm_neighbors`, and `dia_neighbors`. `RunOptions` adds
`frame_batch_size` and `skip_validation`; use
`..RunOptions::default()` when constructing options. The canonical empty codec
record now decodes successfully. Reports use schema version 1. The GUI requires
Rust 1.95+, while the root CLI/library minimum remains 1.85.
`Acquisition` adds `PrmPasef` and `Mixed`; update exhaustive matches. Detection
now validates PRM metadata and classifies multiple nonzero frame types as mixed
instead of choosing DDA or DIA by precedence.
`ActiveGates` adds `prm_per_event`; use its default when constructing partial
values. Reports retain schema version 1 with this additive field.


