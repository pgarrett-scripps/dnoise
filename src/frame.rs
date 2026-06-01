//! Flat per-point view of a timsTOF frame in native integer `(scan, tof, intensity)` space.

/// A frame expanded into three parallel per-point arrays plus its scan count.
///
/// `timsrust::Frame` stores points in CSR form (`scan_offsets` is a row pointer
/// into `tof_indices` / `intensities`). The filter works on flat per-point
/// arrays, so we expand once on load and regroup once before encoding.
#[derive(Debug, Clone)]
pub struct FlatFrame {
    /// Frame `Id` from the `Frames` table (used for the DB update WHERE clause).
    pub frame_id: usize,
    /// Number of ion-mobility scans in this frame.
    pub num_scans: usize,
    /// Per-point scan index (`0..num_scans`).
    pub scan: Vec<u32>,
    /// Per-point TOF index.
    pub tof: Vec<u32>,
    /// Per-point intensity.
    pub intensity: Vec<u32>,
}

impl FlatFrame {
    /// Expand a `timsrust::Frame` into flat per-point arrays.
    pub fn from_frame(frame: &timsrust::Frame) -> Self {
        let num_scans = frame.scan_offsets.len().saturating_sub(1);
        let n = frame.tof_indices.len();
        let mut scan = Vec::with_capacity(n);
        for s in 0..num_scans {
            let count = frame.scan_offsets[s + 1] - frame.scan_offsets[s];
            scan.extend(std::iter::repeat_n(s as u32, count));
        }
        Self {
            frame_id: frame.index,
            num_scans,
            scan,
            tof: frame.tof_indices.clone(),
            intensity: frame.intensities.clone(),
        }
    }

    /// Number of points in the frame.
    pub fn len(&self) -> usize {
        self.tof.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tof.is_empty()
    }

    /// Collect the points selected by `keep` into `(scan, tof, intensity)` triples,
    /// preserving point order.
    pub fn survivors(&self, keep: &[bool]) -> Vec<(u32, u32, u32)> {
        (0..self.len())
            .filter(|&i| keep[i])
            .map(|i| (self.scan[i], self.tof[i], self.intensity[i]))
            .collect()
    }
}
