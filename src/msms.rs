//! ddaPASEF MS/MS denoising (precursor-centric). The algorithm lives in
//! `dnoise_core::msms` (re-exported here); this module adds the reader loop
//! that feeds every MS/MS frame to the keep-set builder and filters the
//! precursors in parallel.

pub use dnoise_core::msms::{
    MsmsKeep, MsmsKeepBuilder, MsmsKeepParts, PrecursorSpectrum, combine_and_filter,
};

use crate::error::{DnoiseError, Result};
use crate::frame::FlatFrame;
use crate::params::{HaloParams, MsmsFilterParams};
use crate::tdf::{FrameMeta, PasefWindow};
use crate::tsr::FrameReader;
use rayon::prelude::*;
use std::collections::HashSet;

/// Build per-precursor keep sets: read each ddaPASEF MS/MS frame once, partition
/// its points to precursors, combine across frames, and filter each precursor's
/// combined spectrum (in parallel).
pub(crate) fn build_msms_keep(
    reader: &FrameReader,
    meta: &[FrameMeta],
    windows: &[PasefWindow],
    params: &MsmsFilterParams,
    halo: Option<&HaloParams>,
    cancel: Option<&std::sync::atomic::AtomicBool>,
) -> Result<MsmsKeep> {
    let mut builder = MsmsKeepBuilder::new(windows);

    // Accumulate a precursor's fragment points across all its frames.
    for (i, m) in meta.iter().enumerate() {
        if cancel.is_some_and(|c| c.load(std::sync::atomic::Ordering::Relaxed)) {
            return Err(DnoiseError::Cancelled);
        }
        if m.num_peaks == 0 || !builder.wants_frame(m.id) {
            continue;
        }
        let frame = reader.get(i).map_err(|e| DnoiseError::FrameRead {
            index: i,
            message: e.to_string(),
        })?;
        builder.add_frame(m.id, &FlatFrame::from_frame(&frame));
    }

    let fp = params.as_filter_params();
    let (spectra, parts) = builder.into_spectra();
    let keep: Vec<HashSet<u64>> = spectra
        .into_par_iter()
        .map(|s| {
            if cancel.is_some_and(|c| c.load(std::sync::atomic::Ordering::Relaxed)) {
                return Err(DnoiseError::Cancelled);
            }
            Ok(s.filter(&fp, halo))
        })
        .collect::<Result<_>>()?;
    Ok(parts.finish(keep))
}
