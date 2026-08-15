//! Tunable parameters for the vertical-IM feature filter (Stage 1 of ALGORITHM.md).

/// Knobs for [`crate::filter`]. Defaults match the tuned PXD070049 benchmark
/// configuration (`benchmark/config/dnoise.toml`).
#[derive(Debug, Clone, Copy)]
pub struct FilterParams {
    /// Column half-width in TOF indices: window spans `[c - w, c + w]`.
    pub mz_half_width: u32,
    /// Minimum number of occupied scans in a kept run (bridged gaps do not count).
    pub min_feature_length: usize,
    /// Max consecutive empty scans tolerated inside a feature (morph-close radius);
    /// bridged scans coalesce neighbouring occupied scans into one run but are not
    /// themselves counted toward `min_feature_length`.
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
            max_internal_gap: 2,
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

/// Knobs for the ddaPASEF MS/MS out-of-window gate: the same scan-interval gate
/// as [`crate::dia_window`], driven by `PasefFrameMsMsInfo` isolation events
/// instead of `DiaFrameMsMsWindows`. On every timsTOF ddaPASEF file examined so
/// far the acquisition writes MS/MS scans only inside the scheduled isolation
/// events, so this gate removes nothing there — it exists as a guarantee, and
/// for acquisitions that behave otherwise.
#[derive(Debug, Clone, Copy, Default)]
pub struct DdaWindowParams {
    /// Scans of leniency added to each side of every isolation event (default 0).
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

/// Knobs for the MS1 selection-polygon gate ([`crate::polygon`]). timsTOF PASEF
/// methods restrict precursor selection to a polygon in the `(m/z, 1/K0)` plane
/// (the "IMS PolygonFilter", stored in `analysis.tdf`). MS1 signal outside it sits
/// where the instrument never schedules fragmentation, so for ddaPASEF it is never
/// a precursor and can be dropped from the survey scans. The polygon itself comes
/// from the data; the pads add physical-unit leniency so a precursor near an edge
/// keeps its isotopic envelope (m/z) and mobility spread (1/K0). Defaults match
/// [`DiaMs1WindowParams`]; set both pads to `0.0` to reproduce the literal
/// polygon.
#[derive(Debug, Clone, Copy)]
pub struct Ms1PolygonParams {
    /// m/z leniency added to each side of the polygon interior, in **Daltons**
    /// (isotopes run to higher m/z, spaced `1/charge` Da).
    pub mz_pad: f64,
    /// Ion-mobility leniency added to each side, in **1/K0**.
    pub im_pad: f64,
}

impl Default for Ms1PolygonParams {
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
    /// Minimum number of occupied scans in a kept run (bridged gaps do not count).
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
            min_feature_length: 3,
            max_internal_gap: 8,
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

/// Physical-unit bounds for the region-of-interest crop ([`crate::crop`]). Every
/// field is optional; an unset bound leaves that side of the axis unconstrained.
/// The crop is a subset of the raw acquisition (not a denoising decision) and, when
/// set, applies to *all* frames (MS1 and MS/MS alike) so a user can carve a smaller
/// `.d` out of a large one. m/z and mobility bounds become integer `(TOF, scan)`
/// ranges via the run calibration; retention time is compared per frame and an
/// out-of-window frame is emitted empty (never deleted), keeping the frame axis and
/// every table that references it valid.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct CropParams {
    /// Lower m/z bound (Da), inclusive. `None` = no lower m/z limit.
    pub mz_min: Option<f64>,
    /// Upper m/z bound (Da), inclusive. `None` = no upper m/z limit.
    pub mz_max: Option<f64>,
    /// Lower ion-mobility bound (`1/K0`), inclusive. `None` = no lower limit.
    pub im_min: Option<f64>,
    /// Upper ion-mobility bound (`1/K0`), inclusive. `None` = no upper limit.
    pub im_max: Option<f64>,
    /// Lower retention-time bound (**minutes**), inclusive. Frames before this are
    /// emitted empty. `None` = no lower RT limit.
    pub rt_min: Option<f64>,
    /// Upper retention-time bound (**minutes**), inclusive. Frames after this are
    /// emitted empty. `None` = no upper RT limit.
    pub rt_max: Option<f64>,
    /// Minimum per-point intensity, inclusive. `None` = no floor.
    pub min_intensity: Option<u32>,
    /// Maximum per-point intensity, inclusive. `None` = no ceiling.
    pub max_intensity: Option<u32>,
}

impl CropParams {
    /// True when no bound is set (the crop would be a no-op).
    pub fn is_empty(&self) -> bool {
        *self == CropParams::default()
    }

    /// True when a retention-time bound is set (frame-level crop).
    pub fn has_rt(&self) -> bool {
        self.rt_min.is_some() || self.rt_max.is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The defaults are load-bearing: they encode the tuned PXD070049 benchmark
    // configuration, and nothing else in the test suite would catch a silent drift.
    #[test]
    fn filter_defaults_match_tuned_benchmark() {
        let p = FilterParams::default();
        assert_eq!(p.mz_half_width, 3);
        assert_eq!(p.min_feature_length, 5);
        assert_eq!(p.max_internal_gap, 2);
        assert_eq!(p.min_window_intensity, 0);
        assert_eq!(p.min_feature_intensity, 0);
        assert_eq!(p.num_iterations, 2);
    }

    #[test]
    fn halo_smooth_watershed_defaults() {
        let h = HaloParams::default();
        assert_eq!(h.peak_fraction, 0.15);
        assert_eq!(h.mz_idx_half_width, 80);
        assert_eq!(h.scan_half_width, 2);

        let s = SmoothParams::default();
        assert_eq!(s.mz_idx_half_width, 2);
        assert_eq!(s.scan_half_width, 3);
        assert_eq!(s.iterations, 1);

        let w = WatershedParams::default();
        assert_eq!(w.box_scan, 10);
        assert_eq!(w.box_mz_idx, 3);
        assert_eq!(w.max_tof_offset, 10);
    }

    #[test]
    fn dia_ms1_defaults_are_physical_pads() {
        let d = DiaMs1WindowParams::default();
        assert_eq!(d.mz_pad, 5.0);
        assert_eq!(d.im_pad, 0.05);
    }

    #[test]
    fn polygon_pad_defaults_match_dia_ms1() {
        // The polygon gate is the ddaPASEF twin of the DIA MS1 gate; its pads
        // default to the same physical leniency.
        let p = Ms1PolygonParams::default();
        let d = DiaMs1WindowParams::default();
        assert_eq!(p.mz_pad, d.mz_pad);
        assert_eq!(p.im_pad, d.im_pad);
    }

    #[test]
    fn msms_defaults_use_a_shorter_feature_length_than_ms1() {
        // Tuned for the short (~25-scan) precursor isolation windows.
        let m = MsmsFilterParams::default();
        assert_eq!(m.min_feature_length, 3);
        assert_eq!(m.max_internal_gap, 8);
        assert_eq!(m.num_iterations, 1);
        assert!(m.min_feature_length < FilterParams::default().min_feature_length);
    }

    #[test]
    fn msms_as_filter_params_preserves_every_knob() {
        let m = MsmsFilterParams {
            mz_half_width: 7,
            min_feature_length: 4,
            max_internal_gap: 9,
            min_window_intensity: 11,
            min_feature_intensity: 13,
            num_iterations: 3,
        };
        let f = m.as_filter_params();
        assert_eq!(f.mz_half_width, m.mz_half_width);
        assert_eq!(f.min_feature_length, m.min_feature_length);
        assert_eq!(f.max_internal_gap, m.max_internal_gap);
        assert_eq!(f.min_window_intensity, m.min_window_intensity);
        assert_eq!(f.min_feature_intensity, m.min_feature_intensity);
        assert_eq!(f.num_iterations, m.num_iterations);
    }

    #[test]
    fn crop_default_is_empty_with_no_rt() {
        let c = CropParams::default();
        assert!(c.is_empty());
        assert!(!c.has_rt());
    }

    #[test]
    fn crop_with_a_non_rt_bound_is_not_empty_but_has_no_rt() {
        let c = CropParams {
            mz_min: Some(400.0),
            ..CropParams::default()
        };
        assert!(!c.is_empty());
        assert!(!c.has_rt());
    }

    #[test]
    fn crop_has_rt_detects_either_bound() {
        assert!(
            CropParams {
                rt_min: Some(1.0),
                ..Default::default()
            }
            .has_rt()
        );
        let hi = CropParams {
            rt_max: Some(9.0),
            ..Default::default()
        };
        assert!(hi.has_rt());
        assert!(!hi.is_empty());
    }

    #[test]
    fn stages_default_is_every_stage_off() {
        let s = Stages::default();
        assert!(!s.filter_all_frames);
        assert_eq!(s.frame_half_width, 0);
        assert!(s.halo.is_none());
        assert!(s.denoise_msms.is_none());
        assert!(s.smooth.is_none());
        assert!(s.watershed.is_none());
        assert!(s.box_centroid.is_none());
        assert!(s.dia_window.is_none());
        assert!(s.dda_window.is_none());
        assert!(!s.dia_per_window);
        assert!(s.dia_ms1.is_none());
        assert!(s.ms1_polygon.is_none());
    }
}

/// Optional pipeline stages layered on top of the core vertical-IM filter
/// ([`FilterParams`]), passed as one value to [`crate::denoise`] and
/// [`crate::denoise_with_progress`] instead of a dozen positional arguments.
///
/// Every field is disabled by default ([`Stages::default`]): the bools are
/// `false`, `frame_half_width` is `0` (off), and the `Option` stages are `None`.
/// Enable a stage by setting its field, borrowing a parameter struct that lives
/// for the call. Each field references the stage documented on its target module.
#[derive(Debug, Default, Clone, Copy)]
pub struct Stages<'a> {
    /// Filter MS/MS frames too. When `false` (default), only MS1 frames are
    /// filtered and MS/MS frames are re-encoded unchanged — the vertical-IM filter
    /// is an MS1 algorithm that strips most MS/MS fragment signal.
    pub filter_all_frames: bool,
    /// Pre-filter MS1 running-average radius (see [`crate::average`]): each MS1
    /// frame's keep/drop decision uses the summed `2*r+1` MS1-frame neighborhood.
    /// `0` (default) reproduces the unsmoothed pipeline. MS/MS frames are never
    /// averaged.
    pub frame_half_width: usize,
    /// Horizontal-halo filter ([`crate::halo`]) after the vertical filter, removing
    /// the weak m/z halo flanking bright peaks. `None` disables it.
    pub halo: Option<&'a HaloParams>,
    /// Denoise MS/MS frames with these knobs instead of passing them through. The
    /// acquisition scheme is auto-detected: ddaPASEF combines each precursor's
    /// fragment scans across frames before filtering (see [`crate::msms`]);
    /// diaPASEF runs the same filter on each whole MS/MS frame. `None` leaves MS/MS
    /// unchanged.
    pub denoise_msms: Option<&'a MsmsFilterParams>,
    /// Box-averaging smoother ([`crate::smooth`]) on each filtered frame's
    /// survivors (after halo, before centroiding) to stabilise watershed seeding.
    /// `None` disables it.
    pub smooth: Option<&'a SmoothParams>,
    /// Watershed centroider ([`crate::watershed`]) as the final stage, collapsing
    /// survivors into intensity-weighted centroids. Mutually exclusive with
    /// `box_centroid`. `None` disables it.
    pub watershed: Option<&'a WatershedParams>,
    /// Greedy small-box centroider ([`crate::box_centroid`]) as the final stage,
    /// tiling streaks into small centroids rather than collapsing them. Mutually
    /// exclusive with `watershed`. `None` disables it.
    pub box_centroid: Option<&'a BoxCentroidParams>,
    /// diaPASEF MS/MS out-of-window gate ([`crate::dia_window`]): drop fragment
    /// points whose mobility scan falls outside every isolation window for their
    /// frame. No effect on ddaPASEF. `None` disables it.
    pub dia_window: Option<&'a DiaWindowParams>,
    /// ddaPASEF MS/MS out-of-window gate: drop fragment points whose mobility
    /// scan falls outside every `PasefFrameMsMsInfo` isolation event for their
    /// frame. Standard timsTOF ddaPASEF acquisitions record no such points, so
    /// this is a guarantee rather than a reduction. No effect on diaPASEF.
    /// `None` disables it.
    pub dda_window: Option<&'a DdaWindowParams>,
    /// diaPASEF per-window MS/MS filtering ([`crate::dia_window::filter_per_window`]):
    /// run the MS/MS filter independently inside each isolation window's scan slice
    /// instead of over the whole frame, so a mobility run cannot fuse across a window
    /// boundary. On by default wherever the MS/MS filter runs; ignored on ddaPASEF.
    pub dia_per_window: bool,
    /// diaPASEF MS1 out-of-window gate ([`crate::dia_ms1`]): drop MS1 points whose
    /// `(m/z, mobility)` falls outside every padded isolation window. No effect on
    /// ddaPASEF. `None` disables it.
    pub dia_ms1: Option<&'a DiaMs1WindowParams>,
    /// ddaPASEF MS1 selection-polygon gate ([`crate::polygon`]): drop MS1 points
    /// outside the run's IMS PolygonFilter (never-selected precursor space).
    /// Auto-detected; skipped on diaPASEF and when the run stores no polygon.
    /// `None` disables it.
    pub ms1_polygon: Option<&'a Ms1PolygonParams>,
}
