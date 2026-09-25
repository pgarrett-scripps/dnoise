//! Bounded temporal evidence for native-point filtering. Observations are
//! indexed once per run ([`NeighborIndexBuilder`]); each frame then sums the
//! points of its compatible neighbors ([`NeighborIndex::keep_mask`]) only to
//! decide which of its own native points survive.
//!
//! **Cross-frame state.** Denoising frame `i` with neighbors enabled reads up to
//! `radius` compatible observations on either side of it, through the caller's
//! frame source. The index itself is read-only after it is built, so frames can
//! be processed in any order and in parallel.

use crate::error::{Error, Result};
use crate::frame::FlatFrame;
use crate::msms::combine_and_filter;
use crate::params::{FilterParams, HaloParams};
use crate::windows::FrameMeta;
use std::collections::HashMap;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

/// Actual temporal evidence for nonempty central events in processed frames.
/// A neighbor is counted only if it supplies points inside the matching event;
/// the central observation is excluded. Counts include repeated uses across events.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize)]
pub struct NeighborUsage {
    /// Nonempty central events evaluated with temporal support enabled.
    pub events: u64,
    /// Evaluated events receiving no points from other observations.
    pub events_without_neighbors: u64,
    /// Total contributing neighbor observations across evaluated events.
    pub neighbors_used: u64,
    /// Largest contributing-neighbor count for one event.
    pub max_neighbors_used: u64,
}

impl NeighborUsage {
    /// Accumulate another frame's usage into this total.
    pub fn add(&mut self, other: Self) {
        self.events += other.events;
        self.events_without_neighbors += other.events_without_neighbors;
        self.neighbors_used += other.neighbors_used;
        self.max_neighbors_used = self.max_neighbors_used.max(other.max_neighbors_used);
    }
}

/// One SQLite-typed metadata value, compared exactly (reals by bit pattern).
/// Frames are only pooled when every compatibility value matches.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Atom {
    /// SQL NULL.
    Null,
    /// INTEGER.
    Integer(i64),
    /// REAL, as `f64::to_bits`.
    Real(u64),
    /// TEXT bytes.
    Text(Vec<u8>),
    /// BLOB bytes.
    Blob(Vec<u8>),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Key {
    target: i64,
    begin: u32,
    end: u32,
    geometry: Arc<[u64]>,
    compatibility: Vec<Atom>,
    segment: usize,
}
#[derive(Clone, Copy, Debug)]
struct Observation {
    frame: usize,
    begin: u32,
    end: u32,
}
#[derive(Debug)]
struct Group {
    radius: usize,
    observations: Vec<Observation>,
}

/// Compatible-observation groups for temporal support. Build with
/// [`NeighborIndexBuilder`].
#[derive(Default, Debug)]
pub struct NeighborIndex {
    groups: Vec<Group>,
    // (group, position), preserving each frame's separate event boundaries.
    frames: Vec<Vec<(usize, usize)>>,
    max_gap: f64,
}

/// Groups observations (a frame's scan interval plus its acquisition identity)
/// into runs of compatible observations, in the order they are added.
#[derive(Debug)]
pub struct NeighborIndexBuilder {
    index: NeighborIndex,
    compatibility: Vec<(Vec<Atom>, usize)>,
    keys: HashMap<Key, usize>,
}

impl NeighborIndexBuilder {
    /// `compatibility[i]` is frame `i`'s compatibility values (e.g. calibration
    /// references, scan mode, polarity, ramp settings) and its calibration
    /// segment number; observations only group with identical values.
    /// `max_rt_gap_seconds` bounds the retention-time distance of a neighbor.
    pub fn new(compatibility: Vec<(Vec<Atom>, usize)>, max_rt_gap_seconds: f64) -> Self {
        Self {
            index: NeighborIndex {
                groups: Vec::new(),
                frames: vec![Vec::new(); compatibility.len()],
                max_gap: max_rt_gap_seconds,
            },
            compatibility,
            keys: HashMap::new(),
        }
    }

    /// Add one observation: scans `[begin, end)` of 0-based frame `frame`, for
    /// `target` (0 when untargeted) with isolation `geometry` (empty for MS1).
    /// It pools with up to `radius` neighbors on either side in its group.
    pub fn add(
        &mut self,
        frame: usize,
        begin: u32,
        end: u32,
        target: i64,
        geometry: Arc<[u64]>,
        radius: usize,
    ) {
        let key = Key {
            target,
            begin,
            end,
            geometry,
            compatibility: self.compatibility[frame].0.clone(),
            segment: self.compatibility[frame].1,
        };
        let index = &mut self.index;
        let group = *self.keys.entry(key).or_insert_with(|| {
            let n = index.groups.len();
            index.groups.push(Group {
                radius,
                observations: Vec::new(),
            });
            n
        });
        let pos = index.groups[group].observations.len();
        index.groups[group]
            .observations
            .push(Observation { frame, begin, end });
        index.frames[frame].push((group, pos));
    }

    /// The finished index.
    pub fn finish(self) -> NeighborIndex {
        self.index
    }
}

impl NeighborIndex {
    /// Scan intervals of frame `frame`'s observations (empty when it has none).
    pub fn intervals(&self, frame: usize) -> Vec<(u32, u32)> {
        self.frames[frame]
            .iter()
            .map(|&(g, p)| {
                let o = self.groups[g].observations[p];
                (o.begin, o.end)
            })
            .collect()
    }

    /// Keep mask for frame `i` (`current`) from its pooled neighborhoods, or
    /// `None` when the frame has no observations. `read_frame(j)` supplies
    /// neighbor frame `j` (0-based); each is read at most once per call.
    /// Neighbors outside `rt_keep`, empty, or beyond the retention-time limit
    /// are skipped.
    // All run-local read/filter inputs are borrowed; no mutable cache is shared
    // across workers, so streaming and parallel writing use identical evidence.
    #[allow(clippy::too_many_arguments)]
    pub fn keep_mask(
        &self,
        read_frame: &dyn Fn(usize) -> Result<FlatFrame>,
        meta: &[FrameMeta],
        i: usize,
        current: &FlatFrame,
        params: &FilterParams,
        halo: Option<&HaloParams>,
        rt_keep: &[bool],
        cancel: Option<&AtomicBool>,
    ) -> Result<Option<(Vec<bool>, NeighborUsage)>> {
        if self.frames[i].is_empty() {
            return Ok(None);
        }
        let mut usage = NeighborUsage::default();
        let mut decoded = HashMap::new();
        let mut mask = vec![false; current.len()];
        for &(g, p) in &self.frames[i] {
            let group = &self.groups[g];
            let event = group.observations[p];
            if !current
                .scan
                .iter()
                .any(|&s| s >= event.begin && s < event.end)
            {
                continue;
            }
            usage.events += 1;
            let mut used = 0;
            let lo = p.saturating_sub(group.radius);
            let hi = p
                .saturating_add(group.radius)
                .min(group.observations.len() - 1);
            let mut points = Vec::new();
            for obs in &group.observations[lo..=hi] {
                if cancel.is_some_and(|c| c.load(Ordering::Relaxed)) {
                    return Err(Error::Cancelled);
                }
                let j = obs.frame;
                if !rt_keep[j]
                    || meta[j].num_peaks == 0
                    || (meta[j].rt - meta[i].rt).abs() > self.max_gap
                {
                    continue;
                }
                let flat = if j == i {
                    current
                } else {
                    if let std::collections::hash_map::Entry::Vacant(e) = decoded.entry(j) {
                        e.insert(read_frame(j)?);
                    }
                    &decoded[&j]
                };
                let before = points.len();
                for k in 0..flat.len() {
                    if flat.scan[k] >= obs.begin && flat.scan[k] < obs.end {
                        points.push((flat.scan[k], flat.tof[k], flat.intensity[k]));
                    }
                }
                if j != i && points.len() > before {
                    used += 1;
                }
            }
            usage.neighbors_used += used;
            usage.max_neighbors_used = usage.max_neighbors_used.max(used);
            usage.events_without_neighbors += u64::from(used == 0);
            let keys = combine_and_filter(
                &points,
                event.begin,
                (event.end - event.begin) as usize,
                params,
                halo,
            );
            for (k, keep) in mask.iter_mut().enumerate() {
                if current.scan[k] >= event.begin && current.scan[k] < event.end {
                    *keep =
                        keys.contains(&(((current.scan[k] as u64) << 32) | current.tof[k] as u64));
                }
            }
        }
        Ok(Some((mask, usage)))
    }
}
