//! Tunable parameters for the vertical-IM feature filter (Stage 1 of ALGORITHM.md).

/// Knobs for [`crate::filter`]. Defaults match the dashboard defaults documented
/// in `ALGORITHM.md` (Parameter reference → Stage 1).
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
            mz_half_width: 2,
            min_feature_length: 5,
            max_internal_gap: 1,
            min_window_intensity: 0,
            min_feature_intensity: 0,
            num_iterations: 1,
        }
    }
}
