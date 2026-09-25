//! dnoise-core — the in-memory denoising stages of
//! [dnoise](https://crates.io/crates/dnoise), for Bruker timsTOF frames.
//!
//! This crate does no I/O and has no native dependencies. A caller reads frames
//! and run metadata itself (for example with `timsrust` and SQLite), then:
//!
//! 1. builds one [`Denoiser`] per run from [`FrameMeta`] and the parameters,
//! 2. attaches any per-run state its stages need (see [`pipeline`]),
//! 3. calls [`Denoiser::process`], [`Denoiser::denoise_ms1`] or
//!    [`Denoiser::denoise_msms`] per frame, in any order, from any thread
//!    (`Denoiser` is `Send + Sync`, checked at compile time).
//!
//! The `dnoise` crate runs exactly this path for its CLI and writes the
//! survivors into a new `.d` folder.
//!
//! # Example
//!
//! ```
//! use dnoise_core::{Acquisition, Denoiser, FilterParams, FlatFrame, FrameMeta, Stages};
//!
//! let meta = vec![FrameMeta { id: 1, num_scans: 700, num_peaks: 8, ms_ms_type: 0, rt: 0.5 }];
//! let params = FilterParams::default();
//! let stages = Stages::default();
//! let denoiser = Denoiser::new(&params, &stages, Acquisition::Ms1Only, meta, false)?;
//!
//! // A vertical streak at TOF 1000 across 8 consecutive mobility scans.
//! let frame = FlatFrame {
//!     frame_id: 1,
//!     num_scans: 700,
//!     scan: (10..18).collect(),
//!     tof: vec![1000; 8],
//!     intensity: vec![100; 8],
//! };
//! let out = denoiser.process(0, &|_| Ok(frame.clone()), None)?;
//! assert_eq!(out.survivors.len(), 8); // the streak survives
//! # Ok::<(), dnoise_core::Error>(())
//! ```
//!
//! The lower-level stages ([`filter::filter_iterated`], [`halo`], the gates,
//! [`smooth`], [`watershed`], [`box_centroid`]) are public too.

#![warn(missing_docs)]

pub mod average;
pub mod box_centroid;
pub mod convert;
pub mod crop;
pub mod dia_ms1;
pub mod dia_window;
pub mod error;
pub mod filter;
pub mod frame;
pub mod halo;
pub mod mobility;
pub mod msms;
pub mod neighbor;
pub mod overlap;
pub mod params;
pub mod pipeline;
pub mod polygon;
pub mod smooth;
pub mod units;
pub mod watershed;
pub mod windows;

pub use error::{Error, Result};
pub use frame::{CsrFrame, FlatFrame};
pub use mobility::{MobilityScale, ScanToMobility};
pub use params::{
    BoxCentroidParams, CropParams, DdaWindowParams, DiaMs1WindowParams, DiaWindowParams,
    FilterParams, HaloParams, Ms1PolygonParams, MsmsFilterParams, NeighborParams, SmoothParams,
    Stages, WatershedParams,
};
pub use pipeline::{Acquisition, DecodedFrame, Denoiser};
pub use windows::{DiaWindows, FrameMeta, PasefWindow};

// Doc-test the README's usage example.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;

// `Denoiser` is shared across threads by embedders; keep it Send + Sync.
const _: () = {
    const fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<Denoiser<'static>>();
};
