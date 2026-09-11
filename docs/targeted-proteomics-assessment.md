# Targeted proteomics support assessment

September 9, 2026 · Investigation of the current local 0.2.0-dev tree.

**Follow-up:** the conservative compatibility stage is now implemented and
tested on one public run. See [current support and evidence](targeted.md).
The investigation below records the behavior before that implementation.

**Recommendation: start with prm-PASEF on timsTOF.** It retains the native
scan/TOF/intensity representation used by dnoise. Reading a targeted run and
preserving its metadata looks achievable; validating fragment denoising is a
separate research step. No production behavior was changed in this investigation.

Conventional PRM and SRM/MRM datasets without resolved ion mobility require a
different filtering model as well as new format readers/writers. They are a
larger extension than adding another PASEF acquisition mode.

## Findings in the current code

| Area | Finding and implication |
|---|---|
| Reading | The pinned timsrust 0.4.2 frame reader accepts acquisition types it labels `Unknown`; it only recognizes 0/8/9 as MS1/DDA/DIA. dnoise reads MS level independently from `Frames.MsMsType`. Basic frame access need not wait for upstream PRM spectrum support. Compression type 2 remains required. |
| Detection | `src/writer.rs` reports only DDA, DIA, or `MS1-only/unknown`. Neither `src/tdf/mod.rs` nor the writer reads PRM event/target tables. |
| Default processing | Nonzero frame types normally retain their decoded MS/MS points when fragment filtering and cropping are off. Metadata tables are copied with the database. This is promising compatibility plumbing, not evidence of validated PRM support. |
| Fragment filtering | Both the file writer and `RunContext::open` fall back to whole-frame MS/MS filtering when DDA events are absent. PRM target boundaries are consequently ignored. |
| Geometry | The MS1 polygon builder excludes DIA geometry, but does not positively require DDA. PRM needs an explicit gate policy. The existing interval representation merges touching windows; that representation must not erase boundaries between distinct PRM events. |

ProteoWizard's reader provides an implementation reference: it recognizes populated
`PrmFrameMsMsInfo`, joins `PrmTargets` through `Target`, and reads the frame,
exclusive scan-end boundary, isolation m/z and width, collision energy, target
monoisotopic m/z, and charge. Table existence alone is insufficient for detection;
inspect populated events and their associated frames. Actual vendor schemas,
frame-type values, and acquisition-software variants still need verification on
real files. [ProteoWizard TimsData.cpp](https://github.com/ProteoWizard/pwiz/blob/master/pwiz_aux/msrc/utility/vendor_api/Bruker/TimsData.cpp)

## Executed compatibility probe

A temporary Rust example reused `tests/common/mod.rs` to generate a four-frame,
type-2 dataset. It changed fragment frame types to 10 (an unknown-type probe,
not independent verification of a vendor enum), added two synthetic PRM events
at scan intervals `[3,7)` and `[7,11)`, and added their target rows. The fragment
frame contained eight equal-intensity points at one TOF index across scans 3–10
plus one isolated point. Halo was disabled and the MS/MS minimum feature length
was set to five scans. The example ran successfully with `cargo run --offline`
and was removed afterward.

| Mode | Result |
|---|---|
| MS1-only filtering | All nine fragment points retained exactly; both event rows retained; output passed full structural/decoding validation. |
| MS/MS filtering enabled | Eight fragment points retained; both event rows retained; output passed full validation. Each target event contained only four occupied scans, demonstrating that the current fallback links the two events into one qualifying feature. |

This establishes a specific dispatch/boundary gap. It does not establish that
discarding either four-scan signal would be scientifically desirable. No real
PRM binary, Bruker SDK round trip, or downstream quantitative comparison was
tested. Structural validation alone cannot establish target-aware filtering.

## Proposed implementation sequence

1. **Compatibility and explicit detection.** Inventory a real run, parse PRM
   events with checked frame/target references and scan bounds, preserve all
   target and scheduling tables, and report PRM in the CLI/GUI. Share dispatch
   between the writer and streaming context. Handle mixed acquisitions per
   frame; reject ambiguous metadata when fragment filtering is requested.
2. **Conservative first support.** Validate MS1-only processing with decoded
   MS/MS equality and metadata equality. Avoid unverified acquisition gates on
   PRM. Report runs with no MS1 frames as having no MS1 denoising opportunity.
   Keep targeted fragment denoising unavailable until its explicit path exists.
3. **Experimental fragment filtering.** Filter individual recorded isolation
   events without merging adjacent targets. Preserve native timestamps and
   survivor intensities. Do not map PRM target IDs directly onto the DDA
   accumulator: it pools all frames for an ID, which could couple decisions
   across a long chromatographic trace. Any temporal support should be bounded,
   target-specific, and separately evaluated.
4. **Quantitative validation.** Compare raw, MS1-only, and experimental outputs
   using identical target lists and extraction settings. Evaluate fragment peak
   areas, chromatographic shape, transition ratios, missingness, replicate CVs,
   and isotope-standard ratios where available. Use a dilution series for
   response linearity and detection/quantification limits. Set acceptance
   tolerances before tuning; retain held-out runs. Measure file reduction only
   alongside those fidelity results.

Allow approximately 3–5 developer days for a first compatibility implementation
once representative files and an independent reader are available. An initial
fragment-denoising evaluation could take 2–4 calendar weeks with suitable
replicates and a dilution series already in hand; these are planning estimates,
not demonstrated performance or a release commitment.

## Public starting dataset

**PXD049405** contains the DIA-to-PRM workflow on a timsTOF Pro. Its publication
reports three technical PRM replicates, making it a candidate for an initial
compatibility and precision comparison. The archive's individual files,
compression types, target tables, and analysis artifacts have not been inspected
or downloaded. It is not established here as a dilution-series benchmark.
[Dataset record](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD049405-1&test=no)
· [Publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC11798701/)

The first implementation prerequisite is inspecting one real PRM `.d` from this
archive or the intended user workflow, including acquisition mode, compression,
target/event schema, calibration segments, and downstream extraction behavior.
