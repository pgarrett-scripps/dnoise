# timsrust 0.4.2 -> 0.6.6 port (branch try/deps-upgrade)

Base: v0.4.0 (9d2054c). Status: `cargo check --workspace --all-targets` and
`cargo clippy --workspace --all-targets -- -D warnings` pass. Nothing has been
run: no tests, no release build, no equivalence run.

## What changed, per file

- `Cargo.toml`: `timsrust = "0.4.2"` replaced by `timsrust-tdf = { version = "0.6.6", default-features = false }`
  (the TDF reader sub-crate only). `rusqlite` 0.32 -> 0.35, the newest that
  resolves: timsrust's `filemanager` pins `rusqlite 0.35` and only one
  `libsqlite3-sys` may be linked.
- `Cargo.lock`: re-resolved offline for the above.
- `src/tsr.rs` (new, `pub mod tsr`): the adapter. It keeps timsrust 0.4.2's
  shapes and formulas so the rest of the crate is unchanged:
  - `ConvertableDomain` with fractional `convert`/`invert`;
  - `Tof2MzConverter` (`mz = (sqrt(mz_min) + slope * tof)^2`, fractional
    inverse) and `Scan2ImConverter` (line from `im_max` at scan 0 to `im_min`
    at `max(NumScans)`, fractional inverse), both 0.4.2's code verbatim;
  - `MetadataReader::new` reading `GlobalMetadata` and `MAX(NumScans)` straight
    from `analysis.tdf` with the same keys, the same otofControl +/-5 m/z
    widening and the same TEXT parsing as 0.4.2 (it also still requires a
    parseable `TimsCompressionType`, as 0.4.2 did);
  - `FrameReader` wrapping `timsrust_tdf::TdfFrameReader`: `get(i)` takes the
    0-based position in sorted `Frames.Id` order and returns a `Frame` with
    plain `index` / `scan_offsets` / `tof_indices: Vec<u32>` /
    `intensities: Vec<u32>`. It rejects `TimsCompressionType != 2` with 0.4.2's
    message and opens via `without_metadata`, so timsrust's per-frame
    `Metadata::new` is never built.
- `src/lib.rs`: `pub mod tsr;`; imports switched to `crate::tsr`.
- `src/writer.rs`, `src/msms.rs`, `src/neighbor.rs`, `src/mobility.rs`,
  `src/frame.rs`: imports / type paths switched from `timsrust::...` to
  `crate::tsr::...`. No logic changes.
- `examples/{bench_filter,check_codec,dump_frame,precursor_area,validate}.rs`:
  imports switched to `dnoise::tsr::...`.
- `tests/pipeline.rs`, `tests/streaming.rs`: `timsrust::readers::FrameReader`
  -> `dnoise::tsr::FrameReader`.
- `tests/targeted.rs`: the `usize`/`u64` SQL params are cast to `i64`
  (needed from rusqlite 0.38; harmless on 0.35).
- No SQL changes were needed for rusqlite 0.35's single-statement rule: every
  `execute`/`prepare`/`query_row` in src, tests and examples is one statement;
  multi-statement SQL already goes through `execute_batch`.
- `dnoise-gui`: no source changes (it only uses the dnoise library).

## Where output could differ from 0.4.0, and why it should not

1. Frame order. 0.4.2 `get(i)` indexed the `Frames` rows in table (rowid)
   order; 0.6.6 indexes by `Frames.Id` and `iter_indices()` is HashMap order.
   The adapter sorts the Ids once and maps position `i` to the i-th Id. `Id` is
   the table's INTEGER PRIMARY KEY (= rowid), so this is the same order, and it
   is the order dnoise's own `tdf::read_frame_meta` (`ORDER BY Id`) uses.
2. m/z conversion. 0.6.6's Tof2Mz inverse truncates to `u32`. The adapter keeps
   0.4.2's converter, so the ppm half-width (`lib.rs`), gate edges and crop
   edges (`writer.rs`) get the same fractional TOF values.
3. Linear 1/K0 (`--linear-mobility`, and the calibrated scale's fallback). 0.6.6
   anchors at `max(NumScans) - 1` and truncates on inverse. The adapter keeps
   0.4.2's line anchored at `max(NumScans)`. The calibrated scale in
   `src/mobility.rs` never used timsrust.
4. Metadata. Same `GlobalMetadata` keys, same TEXT -> number parsing, same
   `MAX(NumScans)` (0.4.2 took the max of the column read as `u32`; SQL `MAX`
   gives the same value on Bruker files, which have no NULL `NumScans`).
5. Frame decoding. timsrust-tdf 0.6.6's type-2 decoder is line-for-line 0.4.2's
   (scan offsets = half the per-scan word counts, TOF = running delta sum - 1,
   zstd `decode_all`). The typed `TofIndex`/`IntensityIndex` store `n + 1` in a
   `NonZeroU32`; `u32::from` gives back `n` exactly. Differences are only in
   failure paths:
   - an empty payload now errors with `EmptyData` (0.4.2: `CorruptFrame`).
     dnoise never reads `NumPeaks == 0` frames through the reader, so this is
     not reached; tests/examples that tolerate an `Err` still do;
   - a payload that decompresses to zero bytes panics in 0.6.6
     (`expect("Blob cannot be empty")`) instead of erroring;
   - a TOF or intensity value of `u32::MAX` panics in 0.6.6 (`n + 1`
     overflow); impossible on real data (TOF < `DigitizerNumSamples`).
6. Frame table parsing. 0.6.6 deserialises `Frames` rows with serde, where
   0.4.2 defaulted unreadable cells to 0; a `.d` with NULLs in `Frames` columns
   could now fail to open. For diaPASEF it also still reads the window groups
   and quadrupole settings (unused by dnoise) at open, as 0.4.2 did.
7. SQLite library. rusqlite 0.35 bundles SQLite 3.49.1 (0.4.0: 3.46.0).
   dnoise copies `analysis.tdf` and updates it in place, and SQLite stamps its
   own version number into the file header (bytes 96-99). So `analysis.tdf`
   will NOT be byte-identical: expect a difference at offset 96-99 only.
   Content (every table) should be identical. `analysis.tdf_bin` does not
   involve SQLite and zstd is unchanged (zstd 0.13.3 / zstd-sys 1.5.6 in both
   locks), so it should be byte-identical.
8. Compression type 1 (older timsTOF files): 0.6.6 can read it, but the adapter
   rejects it as 0.4.2 did, so behaviour is unchanged.

## Equivalence test to run later

Use the same benchmark inputs as the 0.4.0 numbers (one ddaPASEF, one diaPASEF
run from the paper benchmark), same flags as the 0.4.0 runs.

1. Build both: `git -C ~/Repos/d_noise worktree` at `v0.4.0` and this branch,
   `cargo build --release` each (queue via ajs, 2-4 cores).
2. Dry run, compare "points kept": 0.4.0 gives DDA 111,722,871 and
   DIA 755,841,088. The port must print the same.
3. Full denoise with both binaries into separate output folders, then:
   - `cmp old/analysis.tdf_bin new/analysis.tdf_bin` (must be identical);
   - `cmp -l old/analysis.tdf new/analysis.tdf` (only bytes 97-100, 1-based,
     i.e. the SQLite version stamp, may differ); and
     `sqlite3 old/analysis.tdf .dump | diff - <(sqlite3 new/analysis.tdf .dump)`
     (must be empty).
4. Repeat 2-3 with `--linear-mobility` (checks the 0.4.2 linear 1/K0 line).
5. Optional: `cargo test --workspace` on this branch.

## Size and dependencies

Counted from `Cargo.lock` and `cargo tree --offline -e normal` (no build):

| | v0.4.0 | this branch |
|---|---|---|
| packages in Cargo.lock (workspace incl. GUI) | 467 | 574 |
| unique crates, `dnoise` CLI (default features) | 170 | 318 |
| unique crates, `dnoise` library only (`--no-default-features`) | 131 | 285 |
| release binary | 6.8 MB (6,810,304 B, `~/Repos/d_noise/target/release/dnoise` 0.4.0) | not built |

The growth is the cloud I/O stack: `timsrust-tdf` -> `timsrust-core` (feature
`io`, required) -> `filemanager` with default features (`cloud`, `parquet`,
`json`, `sql`) -> `object_store` (aws/azure/gcp/http), `reqwest`, `hyper`,
`rustls`, `tokio`, `arrow` 57, `parquet` 57. `timsrust-core` declares
`filemanager` with `default-features = true`, and every reader crate
(`timsrust-tdf`, `-minitdf`, `-tsf`, `-parquet-spectra`) needs `io`, so no
dependency declaration in dnoise can drop it. Options if that matters: a
`[patch.crates-io]` fork of `timsrust-core` with `default-features = false` on
`filemanager`, an upstream fix, or decoding frames with dnoise's own
`codec::decode` and dropping timsrust (the decoder is ~30 lines and already
mirrors timsrust's).
