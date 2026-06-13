//! dnoise — denoise Bruker timsTOF `.d` folders with the iterative vertical-IM
//! feature filter (Stage 1 of ALGORITHM.md). Reads via `timsrust`, re-encodes
//! surviving points as a type-2 `tdf_bin`, and writes a new `.d` folder.
//!
//! Two API tiers are exposed:
//!
//! - **High level** — [`denoise`] / [`denoise_with_progress`] take an input `.d`
//!   and write a filtered output `.d`, parameterized by [`FilterParams`].
//! - **Low level** — [`FlatFrame`] plus [`filter::filter_once`] /
//!   [`filter::filter_iterated`], [`average::running_average`], and the type-2
//!   [`codec`] operate on in-memory frames for callers that bring their own I/O.
//!
//! # Example
//!
//! Run the vertical-IM filter over an in-memory frame:
//!
//! ```
//! use dnoise::{FilterParams, FlatFrame, filter_iterated};
//!
//! // A vertical streak at TOF 1000 across 8 consecutive mobility scans.
//! let frame = FlatFrame {
//!     frame_id: 1,
//!     num_scans: 700,
//!     scan: (10..18).collect(),
//!     tof: vec![1000; 8],
//!     intensity: vec![100; 8],
//! };
//! let keep = filter_iterated(&frame, &FilterParams::default());
//! assert!(keep.iter().all(|&k| k)); // the streak survives
//! ```

#![warn(missing_docs)]

pub mod average;
pub mod box_centroid;
pub mod codec;
pub mod dia_ms1;
pub mod dia_window;
pub mod error;
pub mod filter;
pub mod frame;
pub mod halo;
pub mod msms;
pub mod params;
pub mod smooth;
pub mod watershed;
pub mod writer;

// SQLite plumbing for the high-level pipeline; not part of the public API.
mod tdf;

// High-level pipeline.
pub use error::{DecodeError, DnoiseError, Result};
pub use params::{
    BoxCentroidParams, DiaMs1WindowParams, DiaWindowParams, FilterParams, HaloParams,
    MsmsFilterParams, SmoothParams, WatershedParams,
};
pub use writer::{DenoiseStats, Progress, denoise, denoise_with_progress};

// Low-level building blocks.
pub use average::running_average;
pub use box_centroid::box_centroid;
pub use dia_ms1::{DiaMs1Gate, TofScanBox};
pub use dia_window::{filter_per_window, in_window_mask};
pub use filter::{filter_iterated, filter_once};
pub use frame::FlatFrame;
pub use halo::horizontal_halo_keep_mask;
pub use msms::combine_and_filter;
pub use smooth::box_average;
pub use watershed::watershed_centroid;
