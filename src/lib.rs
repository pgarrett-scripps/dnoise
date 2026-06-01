//! dnoise — denoise Bruker timsTOF `.d` folders with the iterative vertical-IM
//! feature filter (Stage 1 of ALGORITHM.md). Reads via `timsrust`, re-encodes
//! surviving points as a type-2 `tdf_bin`, and writes a new `.d` folder.

pub mod average;
pub mod filter;
pub mod frame;
pub mod params;
pub mod tdf;
pub mod writer;

pub use params::FilterParams;
pub use writer::{DenoiseStats, denoise};
