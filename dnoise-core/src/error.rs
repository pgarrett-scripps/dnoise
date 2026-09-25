//! The error type of the pure stages and the per-frame API.

use thiserror::Error;

/// Convenience alias for results produced by this crate.
pub type Result<T, E = Error> = std::result::Result<T, E>;

/// Anything that can go wrong while building the run state or denoising a frame.
/// The messages match the `dnoise` crate's `DnoiseError` variants of the same
/// name, which this converts into one to one.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// Invalid input, configuration or metadata.
    #[error("{0}")]
    InvalidInput(String),

    /// Reading or decoding a supporting frame failed (reported by the caller's
    /// frame source).
    #[error("reading frame {index}: {message}")]
    FrameRead {
        /// Zero-based frame index that failed.
        index: usize,
        /// Upstream reader message.
        message: String,
    },

    /// A cancellation token was set before the work finished.
    #[error("cancelled")]
    Cancelled,
}
