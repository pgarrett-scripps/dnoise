//! Tunable parameters for the vertical-IM feature filter (Stage 1 of ALGORITHM.md).

/// Knobs for [`crate::filter`]. Defaults match the tuned PXD070049 benchmark
/// configuration (`benchmark/config/dnoise.toml`).
#[derive(Debug, Clone, Copy)]
pub struct FilterParams {
    /// Column half-width in TOF indices: window spans `[c - w, c + w]`.
    pub mz_half_width: u32,
    /// Minimum total span (gap-inclusive, in scans) of a kept run.
    pub min_feature_length: usize,
    /// Max consecutive empty scans tolerated inside a feature (morph-close radius).
    pub max_internal_gap: usize,
    /// Per-scan summed-intensity floor for marking a scan occupied.
    pub min_window_intensity: u64,
    /// Total summed intensity (window × span) required for a run to be kept.
    pub min_feature_intensity: u64,
    /// How many times to re-apply the filter to its own survivors.
    pub num_iterations: usize,
}

impl Default for FilterParams {
    fn default() -> Self {
        Self {
            mz_half_width: 3,
            min_feature_length: 5,
            max_internal_gap: 1,
            min_window_intensity: 0,
            min_feature_intensity: 0,
            num_iterations: 2,
        }
    }
}

/// Knobs for the horizontal-halo filter ([`crate::halo`]), which removes the
/// weak m/z halo flanking bright ions. Each point is compared to the maximum
/// intensity of its surrounding box *excluding its own TOF column*, and dropped
/// if below `peak_fraction` of that reference. Operates in integer
/// `(scan, TOF index)` space. Defaults match the tuned PXD070049 benchmark
/// configuration (`benchmark/config/dnoise.toml`).
#[derive(Debug, Clone, Copy)]
pub struct HaloParams {
    /// Drop a peak below this fraction of the off-column box-max reference.
    pub peak_fraction: f64,
    /// Half-width of the reference box along TOF index.
    pub mz_idx_half_width: u32,
    /// Half-width of the reference box along ion-mobility scan.
    pub scan_half_width: usize,
}

impl Default for HaloParams {
    fn default() -> Self {
        Self {
            peak_fraction: 0.15,
            mz_idx_half_width: 80,
            scan_half_width: 2,
        }
    }
}

/// Knobs for the box-averaging smoother ([`crate::smooth`]). Rewrites each
/// surviving point's intensity with the mean intensity of all points inside a
/// `(±scan_half_width, ±mz_idx_half_width)` box (inclusive of the point itself),
/// in integer `(scan, TOF index)` space — coordinates are unchanged. Runs after
/// the halo filter and before the watershed centroider to stabilise seeding:
/// noise-driven local intensity spikes otherwise make watershed split one ion
/// into several centroids.
#[derive(Debug, Clone, Copy)]
pub struct SmoothParams {
    /// Half-width of the averaging box along TOF index.
    pub mz_idx_half_width: u32,
    /// Half-width of the averaging box along ion-mobility scan.
    pub scan_half_width: usize,
    /// How many times to re-apply the smoother to its own output.
    pub iterations: usize,
}

impl Default for SmoothParams {
    fn default() -> Self {
        Self {
            mz_idx_half_width: 2,
            scan_half_width: 3,
            iterations: 1,
        }
    }
}

/// Knobs for the watershed centroider ([`crate::watershed`]), an optional final
/// stage that collapses each watershed group of raw points into a single
/// intensity-weighted `(scan, TOF index)` centroid. Unlike the vertical and halo
/// filters (which select a subset of points), this is a *lossy* reduction that
/// typically shrinks the surviving point count to a small fraction. Operates in
/// integer `(scan, TOF index)` space. Defaults match `koth_rust`'s watershed.
#[derive(Debug, Clone, Copy)]
pub struct WatershedParams {
    /// Nearest-neighbour reach along the ion-mobility scan axis.
    pub box_scan: u32,
    /// Nearest-neighbour reach along the TOF-index axis.
    pub box_mz_idx: u32,
    /// Minimum intensity for a point to open (seed) a new group.
    pub min_seed_intensity: u64,
    /// Drop groups whose summed intensity is below this.
    pub min_centroid_total: u64,
    /// Hard cap on how far a follower may sit from its group's seed in TOF
    /// indices — prevents a long follower chain from creeping past the peak edge.
    pub max_tof_offset: u32,
}

impl Default for WatershedParams {
    fn default() -> Self {
        Self {
            box_scan: 10,
            box_mz_idx: 3,
            min_seed_intensity: 0,
            min_centroid_total: 0,
            max_tof_offset: 10,
        }
    }
}

/// Knobs for the greedy small-box centroider ([`crate::box_centroid`]), an
/// optional final reduction stage. Each box `(±scan_half_width, ±mz_idx_half_width)`
/// of points is merged into one intensity-weighted `(scan, TOF index)` centroid;
/// unlike the watershed centroider this is non-transitive, so a long mobility
/// streak is tiled into several small centroids (preserving its profile) rather
/// than collapsed to a point. Boxes with summed intensity below
/// `min_centroid_total` are dropped (the optional denoising floor; 0 conserves
/// total intensity exactly).
#[derive(Debug, Clone, Copy)]
pub struct BoxCentroidParams {
    /// Box half-width along TOF index (m/z); keep tight to preserve m/z precision.
    pub mz_idx_half_width: u32,
    /// Box half-width along ion-mobility scan.
    pub scan_half_width: u32,
    /// Drop a box whose summed intensity is below this (0 = keep all).
    pub min_centroid_total: u64,
}

impl Default for BoxCentroidParams {
    fn default() -> Self {
        Self {
            mz_idx_half_width: 2,
            scan_half_width: 2,
            min_centroid_total: 0,
        }
    }
}

/// Knobs for the diaPASEF isolation-window MS/MS filter ([`crate::dia_window`]).
/// The window scan intervals come from the data (`DiaFrameMsMsWindows`); the only
/// tunable is `scan_pad`, which symmetrically widens each window to tolerate
/// signal a few mobility scans past an isolation edge before it is gated out.
#[derive(Debug, Clone, Copy, Default)]
pub struct DiaWindowParams {
    /// Scans of leniency added to each side of every isolation window (default 0).
    pub scan_pad: u32,
}

/// Knobs for the diaPASEF **MS1** out-of-window gate ([`crate::dia_ms1`]). DIA
/// tiles the precursor space into isolation windows, each covering an m/z band
/// over a mobility-scan interval (`DiaFrameMsMsWindows`). An MS1 peak that falls
/// in no window's `(m/z, mobility)` region is a precursor that is never isolated,
/// so it is dropped. Each window is padded — in **physical** units, since the
/// windows are defined in m/z and the calibration is needed anyway — so a
/// precursor near a window edge keeps its full isotopic envelope (isotopes run to
/// higher m/z) and mobility spread.
#[derive(Debug, Clone, Copy)]
pub struct DiaMs1WindowParams {
    /// m/z leniency added to each side of every window, in **Daltons**. Maps
    /// directly to isotopes (spaced `1/charge` Da), uniformly across the m/z range.
    pub mz_pad: f64,
    /// Ion-mobility leniency added to each side of every window, in **1/K0**.
    pub im_pad: f64,
}

impl Default for DiaMs1WindowParams {
    fn default() -> Self {
        Self {
            mz_pad: 5.0,
            im_pad: 0.05,
        }
    }
}

/// Vertical-filter knobs for the ddaPASEF MS/MS path ([`crate::msms`]). Mirrors
/// [`FilterParams`] but with defaults tuned for the short (~25-scan) precursor
/// isolation windows: a smaller `min_feature_length`.
#[derive(Debug, Clone, Copy)]
pub struct MsmsFilterParams {
    /// Column half-width in TOF indices.
    pub mz_half_width: u32,
    /// Minimum total span (scans) of a kept run.
    pub min_feature_length: usize,
    /// Max consecutive empty scans tolerated inside a feature.
    pub max_internal_gap: usize,
    /// Per-scan summed-intensity floor for occupancy.
    pub min_window_intensity: u64,
    /// Total summed-intensity floor for a kept run.
    pub min_feature_intensity: u64,
    /// Filter passes over the combined spectrum's survivors.
    pub num_iterations: usize,
}

impl Default for MsmsFilterParams {
    fn default() -> Self {
        Self {
            mz_half_width: 3,
            min_feature_length: 2,
            max_internal_gap: 5,
            min_window_intensity: 0,
            min_feature_intensity: 0,
            num_iterations: 1,
        }
    }
}

impl MsmsFilterParams {
    /// View these knobs as a [`FilterParams`] for the vertical filter.
    pub fn as_filter_params(&self) -> FilterParams {
        FilterParams {
            mz_half_width: self.mz_half_width,
            min_feature_length: self.min_feature_length,
            max_internal_gap: self.max_internal_gap,
            min_window_intensity: self.min_window_intensity,
            min_feature_intensity: self.min_feature_intensity,
            num_iterations: self.num_iterations,
        }
    }
}
