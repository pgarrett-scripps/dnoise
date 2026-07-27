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
#[cfg(feature = "config")]
pub mod config;
pub mod crop;
pub mod dia_ms1;
pub mod dia_window;
pub mod error;
pub mod filter;
pub mod frame;
pub mod halo;
pub mod msms;
pub mod params;
pub mod polygon;
pub mod smooth;
pub mod watershed;
pub mod writer;

// SQLite plumbing for the high-level pipeline; not part of the public API.
mod tdf;

// High-level pipeline.
pub use error::{DecodeError, DnoiseError, Result};
pub use params::{
    BoxCentroidParams, CropParams, DiaMs1WindowParams, DiaWindowParams, FilterParams, HaloParams,
    Ms1PolygonParams, MsmsFilterParams, SmoothParams, Stages, WatershedParams,
};
pub use writer::{
    DecodedFrame, DenoiseStats, FrameCtx, Progress, RunOptions, SampleSpec, denoise,
    denoise_with_options, denoise_with_progress, process_frame_decoded,
};

// Low-level building blocks.
pub use crop::CropGate;

use std::path::Path;
use timsrust::converters::ConvertableDomain;
use timsrust::readers::MetadataReader;

/// Acquisition scheme of a `.d` run, detected from the `Frames.MsMsType` column.
/// Drives the `--preset auto` gate selection (see the CLI): ddaPASEF wants the MS1
/// selection-polygon gate, diaPASEF the isolation-window gates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Acquisition {
    /// ddaPASEF (`MsMsType` 8 present): data-dependent PASEF.
    DdaPasef,
    /// diaPASEF (`MsMsType` 9 present): data-independent PASEF.
    DiaPasef,
    /// Only MS1 frames — no MS/MS in the run.
    Ms1Only,
    /// MS/MS frames present but of an unrecognised `MsMsType`.
    Unknown,
}

/// Detect the acquisition scheme of an input `.d` folder by inspecting the frame
/// table (`MsMsType`): 9 anywhere = diaPASEF, else 8 = ddaPASEF, else MS1-only when
/// every frame is MS1, else unknown. A cheap SQLite read, no frame decoding.
pub fn detect_acquisition(input: &Path) -> Result<Acquisition> {
    let meta = tdf::read_frame_meta(&input.join("analysis.tdf"))?;
    Ok(if meta.iter().any(|m| m.ms_ms_type == 9) {
        Acquisition::DiaPasef
    } else if meta.iter().any(|m| m.ms_ms_type == 8) {
        Acquisition::DdaPasef
    } else if meta.iter().all(|m| m.is_ms1()) {
        Acquisition::Ms1Only
    } else {
        Acquisition::Unknown
    })
}

/// Convert a mass tolerance in ppm to a TOF-index half-width for the vertical
/// filter's `mz_half_width`, evaluated at a reference m/z. The timsTOF calibration
/// maps m/z to a fractional TOF index nonlinearly, so a fixed ppm corresponds to a
/// different index width across the mass range; the vertical filter uses one
/// constant index window, so we evaluate the width at `ref_mz` (defaulting to the
/// midpoint of the acquired m/z range) and round up to at least 1.
///
/// `ref_mz`: `Some(mz)` to fix the reference, or `None` to use the acquisition
/// midpoint (falling back to 800.0 if the range is unavailable).
pub fn tof_half_width_for_ppm(input: &Path, ppm: f64, ref_mz: Option<f64>) -> Result<u32> {
    let in_tdf = input.join("analysis.tdf");
    let md = MetadataReader::new(&in_tdf).map_err(|e| DnoiseError::Metadata(e.to_string()))?;
    let mz = match ref_mz {
        Some(m) => m,
        None => match tdf::read_mz_acq_range(&in_tdf)? {
            Some((lo, hi)) => 0.5 * (lo + hi),
            None => 800.0,
        },
    };
    // Half the full ppm span, converted to TOF indices at the reference m/z.
    let t_center = md.mz_converter.invert(mz);
    let t_edge = md.mz_converter.invert(mz * (1.0 + ppm * 1e-6));
    let half = (t_edge - t_center).abs().round() as u32;
    Ok(half.max(1))
}

// Low-level building blocks.
pub use average::running_average;
pub use box_centroid::box_centroid;
pub use dia_ms1::{DiaMs1Gate, TofScanBox};
pub use dia_window::{filter_per_window, in_window_mask};
pub use filter::{filter_iterated, filter_once};
pub use frame::FlatFrame;
pub use halo::horizontal_halo_keep_mask;
pub use msms::combine_and_filter;
pub use polygon::PolygonGate;
pub use smooth::box_average;
pub use watershed::watershed_centroid;
