//! Flat per-point view of a timsTOF frame in native integer `(scan, tof, intensity)` space.

/// A decoded frame in CSR layout: `scan_offsets` is a row pointer, so scan `s`
/// holds points `scan_offsets[s]..scan_offsets[s + 1]` of `tof_indices` /
/// `intensities`. This is the layout timsrust decodes to; implement it for your
/// reader's frame type to use [`FlatFrame::from_frame`], or call
/// [`FlatFrame::from_csr`] with the slices directly.
pub trait CsrFrame {
    /// Frame identifier (`Frames.Id`).
    fn frame_id(&self) -> usize;
    /// CSR row pointer with `num_scans + 1` entries.
    fn scan_offsets(&self) -> &[usize];
    /// Per-point TOF index.
    fn tof_indices(&self) -> &[u32];
    /// Per-point raw intensity.
    fn intensities(&self) -> &[u32];
}

/// A frame expanded into three parallel per-point arrays plus its scan count.
///
/// Readers decode frames in CSR form ([`CsrFrame`]: `scan_offsets` is a row
/// pointer into `tof_indices` / `intensities`). The filter works on flat
/// per-point arrays, so we expand once on load and regroup once before encoding.
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
    /// Expand a CSR frame ([`CsrFrame`]) into flat per-point arrays.
    pub fn from_frame<F: CsrFrame + ?Sized>(frame: &F) -> Self {
        Self::from_csr(
            frame.frame_id(),
            frame.scan_offsets(),
            frame.tof_indices(),
            frame.intensities(),
        )
    }

    /// Expand CSR slices into flat per-point arrays. `scan_offsets` has
    /// `num_scans + 1` entries; `tof_indices` and `intensities` are parallel.
    pub fn from_csr(
        frame_id: usize,
        scan_offsets: &[usize],
        tof_indices: &[u32],
        intensities: &[u32],
    ) -> Self {
        let num_scans = scan_offsets.len().saturating_sub(1);
        let n = tof_indices.len();
        let mut scan = Vec::with_capacity(n);
        for s in 0..num_scans {
            let count = scan_offsets[s + 1] - scan_offsets[s];
            scan.extend(std::iter::repeat_n(s as u32, count));
        }
        Self {
            frame_id,
            num_scans,
            scan,
            tof: tof_indices.to_vec(),
            intensity: intensities.to_vec(),
        }
    }

    /// Number of points in the frame.
    pub fn len(&self) -> usize {
        self.tof.len()
    }

    /// True when the frame contains no points.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct Csr {
        index: usize,
        scan_offsets: Vec<usize>,
        tof_indices: Vec<u32>,
        intensities: Vec<u32>,
    }

    impl CsrFrame for Csr {
        fn frame_id(&self) -> usize {
            self.index
        }
        fn scan_offsets(&self) -> &[usize] {
            &self.scan_offsets
        }
        fn tof_indices(&self) -> &[u32] {
            &self.tof_indices
        }
        fn intensities(&self) -> &[u32] {
            &self.intensities
        }
    }

    fn frame(scan: Vec<u32>, tof: Vec<u32>, intensity: Vec<u32>) -> FlatFrame {
        FlatFrame {
            frame_id: 0,
            num_scans: 4,
            scan,
            tof,
            intensity,
        }
    }

    #[test]
    fn len_and_is_empty() {
        assert_eq!(frame(vec![], vec![], vec![]).len(), 0);
        assert!(frame(vec![], vec![], vec![]).is_empty());

        let g = frame(vec![0, 1], vec![10, 20], vec![1, 2]);
        assert_eq!(g.len(), 2);
        assert!(!g.is_empty());
    }

    #[test]
    fn survivors_selects_kept_points_in_order() {
        let f = frame(vec![0, 1, 2], vec![100, 200, 300], vec![5, 6, 7]);
        assert_eq!(
            f.survivors(&[true, false, true]),
            vec![(0, 100, 5), (2, 300, 7)]
        );
    }

    #[test]
    fn survivors_empty_when_nothing_kept() {
        let f = frame(vec![0], vec![1], vec![1]);
        assert!(f.survivors(&[false]).is_empty());
    }

    #[test]
    fn from_frame_expands_csr_offsets_into_per_point_scans() {
        // `scan_offsets` is a CSR row pointer over 3 scans: scan 0 has two points,
        // scan 1 has none, scan 2 has one.
        let src = Csr {
            scan_offsets: vec![0, 2, 2, 3],
            tof_indices: vec![10, 11, 12],
            intensities: vec![100, 101, 102],
            index: 7,
        };
        let flat = FlatFrame::from_frame(&src);
        assert_eq!(flat.frame_id, 7);
        assert_eq!(flat.num_scans, 3);
        assert_eq!(flat.scan, vec![0, 0, 2]);
        assert_eq!(flat.tof, vec![10, 11, 12]);
        assert_eq!(flat.intensity, vec![100, 101, 102]);
        assert_eq!(flat.len(), 3);
    }

    #[test]
    fn from_frame_handles_an_empty_frame() {
        let src = Csr {
            scan_offsets: vec![0, 0, 0],
            ..Default::default()
        };
        let flat = FlatFrame::from_frame(&src);
        assert_eq!(flat.num_scans, 2);
        assert!(flat.is_empty());
    }
}
