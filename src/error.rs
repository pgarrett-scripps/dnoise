//! Error types for the `dnoise` library.
//!
//! The library returns the typed [`DnoiseError`]; binaries can convert it into
//! `anyhow`/`eyre`/etc. via the standard [`std::error::Error`] impl. The low-level
//! type-2 decoder returns the narrower [`DecodeError`], which folds into
//! [`DnoiseError`] when used through the high-level pipeline.

use std::path::PathBuf;
use thiserror::Error;

/// Convenience alias for results produced by this crate.
pub type Result<T, E = DnoiseError> = std::result::Result<T, E>;

/// Anything that can go wrong while denoising a `.d` folder.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum DnoiseError {
    /// The input path is not a Bruker `.d` folder (missing `analysis.tdf` /
    /// `analysis.tdf_bin`).
    #[error("{0} is not a Bruker .d folder (missing analysis.tdf / analysis.tdf_bin)")]
    NotADotD(PathBuf),

    /// The output folder already exists and overwriting was not requested.
    #[error("output {0} already exists (pass force to overwrite)")]
    OutputExists(PathBuf),

    /// Opening the timsTOF frame data failed.
    #[error("opening frame data: {0}")]
    OpenFrames(String),

    /// Reading the timsTOF metadata / calibration failed (needed by the diaPASEF
    /// MS1 window gate to convert isolation m/z to TOF indices).
    #[error("reading metadata: {0}")]
    Metadata(String),

    /// Reading or decoding a single frame failed.
    #[error("reading frame {index}: {message}")]
    FrameRead {
        /// Zero-based frame index that failed.
        index: usize,
        /// Upstream reader message.
        message: String,
    },

    /// Decoding a type-2 frame record failed.
    #[error(transparent)]
    Decode(#[from] DecodeError),

    /// An underlying I/O error (copying the `.d`, writing `analysis.tdf_bin`, …).
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// An error from the `analysis.tdf` SQLite database.
    #[error("database error: {0}")]
    Sqlite(#[from] rusqlite::Error),

    /// The run was cancelled via the [`crate::RunOptions::cancel`] token before it
    /// finished. Any partial output is incomplete and should be discarded.
    #[error("cancelled")]
    Cancelled,
}

/// Failure decoding a type-2 `tdf_bin` frame record (see [`crate::codec`]).
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum DecodeError {
    /// Record is shorter than the fixed 8-byte header.
    #[error("record shorter than header")]
    ShortRecord,

    /// The record's declared `total_byte_count` is inconsistent with its length.
    #[error("invalid total_byte_count {0}")]
    InvalidByteCount(usize),

    /// The decompressed payload is not a whole number of `u32` values.
    #[error("decompressed blob not u32-aligned")]
    Misaligned,

    /// The decoded `scan_count` exceeds the decompressed payload length.
    #[error("scan_count {scan_count} exceeds blob length {len}")]
    ScanCountOverflow {
        /// Decoded scan count.
        scan_count: usize,
        /// Number of `u32` values actually present.
        len: usize,
    },

    /// zstd decompression of the payload failed.
    #[error("zstd decompression failed: {0}")]
    Zstd(#[source] std::io::Error),
}
