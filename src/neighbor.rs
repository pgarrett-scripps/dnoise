//! Bounded temporal evidence for native-point filtering. Metadata is indexed once;
//! each worker decodes only its selected neighbors, never an entire run of spectra.
use crate::msms::combine_and_filter;
use crate::provenance::NeighborUsage;
use crate::tdf::FrameMeta;
use crate::tsr::FrameReader;
use crate::{Acquisition, DnoiseError, FilterParams, FlatFrame, HaloParams, Result, Stages};
use rusqlite::{Connection, OpenFlags, types::ValueRef};
use std::collections::HashMap;
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Atom {
    Null,
    Integer(i64),
    Real(u64),
    Text(Vec<u8>),
    Blob(Vec<u8>),
}
impl From<ValueRef<'_>> for Atom {
    fn from(v: ValueRef<'_>) -> Self {
        match v {
            ValueRef::Null => Self::Null,
            ValueRef::Integer(x) => Self::Integer(x),
            ValueRef::Real(x) => Self::Real(x.to_bits()),
            ValueRef::Text(x) => Self::Text(x.to_vec()),
            ValueRef::Blob(x) => Self::Blob(x.to_vec()),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Key {
    target: i64,
    begin: u32,
    end: u32,
    geometry: std::sync::Arc<[u64]>,
    compatibility: Vec<Atom>,
    segment: usize,
}
#[derive(Clone, Copy)]
struct Observation {
    frame: usize,
    begin: u32,
    end: u32,
}
struct Group {
    radius: usize,
    observations: Vec<Observation>,
}
#[derive(Default)]
pub(crate) struct NeighborIndex {
    groups: Vec<Group>,
    // (group, position), preserving each frame's separate event boundaries.
    frames: Vec<Vec<(usize, usize)>>,
    max_gap: f64,
}
fn invalid(message: impl Into<String>) -> DnoiseError {
    DnoiseError::InvalidInput(format!("neighbor support: {}", message.into()))
}

impl NeighborIndex {
    pub(crate) fn build(
        path: &Path,
        meta: &[FrameMeta],
        scheme: Acquisition,
        stages: &Stages,
    ) -> Result<Option<Self>> {
        let ms1 = stages.frame_half_width;
        let msms = match scheme {
            Acquisition::PrmPasef => stages.neighbors.prm_radius,
            Acquisition::DiaPasef => stages.neighbors.dia_radius,
            _ => 0,
        };
        if ms1 == 0 && msms == 0 {
            return Ok(None);
        }
        if msms > 0 && stages.denoise_msms.is_none() && !stages.filter_all_frames {
            return Err(invalid(
                "PRM/DIA radius requires MS/MS denoising to be enabled",
            ));
        }
        let db = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
        let columns: Vec<String> = db
            .prepare("PRAGMA table_info(Frames)")?
            .query_map([], |r| r.get(1))?
            .collect::<std::result::Result<_, _>>()?;
        // Pooling raw coordinates requires known calibration references.
        for column in ["MzCalibration", "TimsCalibration"] {
            if !columns.iter().any(|c| c == column) {
                return Err(invalid(format!("Frames.{column} is required")));
            }
        }
        let optional = ["ScanMode", "Polarity", "AccumulationTime", "RampTime"];
        let fields = optional
            .iter()
            .map(|c| {
                if columns.iter().any(|x| x == c) {
                    *c
                } else {
                    "NULL"
                }
            })
            .collect::<Vec<_>>()
            .join(",");
        let measurement_mode: bool = db.query_row("SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type='table' AND name='PrmFrameMeasurementMode')", [], |r| r.get(0))?;
        let mode = if measurement_mode {
            "(SELECT MeasurementModeId FROM PrmFrameMeasurementMode p WHERE p.Frame=f.Id)"
        } else {
            "NULL"
        };
        let mut stmt = db.prepare(&format!("SELECT Id,MzCalibration,TimsCalibration,NumScans,MsMsType,{fields},{mode} FROM Frames f ORDER BY Id"))?;
        let count = stmt.column_count();
        let mut rows = stmt.query([])?;
        let mut compatibility = Vec::with_capacity(meta.len());
        let mut segment = 0;
        let mut previous_cal = None;
        let mut previous_rt = f64::NEG_INFINITY;
        for frame in meta {
            if !frame.rt.is_finite() || frame.rt < previous_rt {
                return Err(invalid(
                    "frame retention times must be finite and nondecreasing",
                ));
            }
            previous_rt = frame.rt;
            let row = rows
                .next()?
                .ok_or_else(|| invalid("frame metadata changed"))?;
            if row.get::<_, usize>(0)? != frame.id {
                return Err(invalid("frame metadata changed"));
            }
            let atoms = (1..count)
                .map(|c| row.get_ref(c).map(Atom::from))
                .collect::<std::result::Result<Vec<_>, _>>()?;
            if matches!(atoms[0], Atom::Null) || matches!(atoms[1], Atom::Null) {
                return Err(invalid("null calibration reference"));
            }
            let cal = (atoms[0].clone(), atoms[1].clone());
            if previous_cal.as_ref().is_some_and(|p| p != &cal) {
                segment += 1;
            }
            previous_cal = Some(cal);
            compatibility.push((atoms, segment));
        }
        let mut index = Self {
            groups: Vec::new(),
            frames: vec![Vec::new(); meta.len()],
            max_gap: stages.neighbors.max_rt_gap_seconds,
        };
        let mut groups = HashMap::<Key, usize>::new();
        let mut add = |frame: usize,
                       begin: u32,
                       end: u32,
                       target: i64,
                       geometry: std::sync::Arc<[u64]>,
                       radius: usize| {
            let key = Key {
                target,
                begin,
                end,
                geometry,
                compatibility: compatibility[frame].0.clone(),
                segment: compatibility[frame].1,
            };
            let group = *groups.entry(key).or_insert_with(|| {
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
        };
        if ms1 > 0 {
            for (i, f) in meta.iter().enumerate().filter(|(_, f)| f.is_ms1()) {
                let end =
                    u32::try_from(f.num_scans).map_err(|_| invalid("scan count exceeds u32"))?;
                add(i, 0, end, 0, Vec::new().into(), ms1);
            }
        }
        if msms > 0 && scheme == Acquisition::DiaPasef {
            let regions = crate::tdf::dia::read(path)?;
            for (i, f) in meta.iter().enumerate() {
                if let Some(group) = regions.frames.get(&f.id) {
                    for r in &regions.groups[group] {
                        add(i, r.begin, r.end, 0, r.signature.clone(), msms);
                    }
                }
            }
        }
        if msms > 0 && scheme == Acquisition::PrmPasef {
            let sql = "SELECT Frame,ScanNumBegin,ScanNumEnd,IsolationMz,IsolationWidth,CollisionEnergy,Target FROM PrmFrameMsMsInfo ORDER BY Frame,ScanNumBegin";
            let expected = 10;
            let by_id: HashMap<_, _> = meta.iter().enumerate().map(|(i, f)| (f.id, i)).collect();
            let mut stmt = db.prepare(sql)?;
            let mut rows = stmt.query([])?;
            let mut ends = HashMap::new();
            while let Some(row) = rows.next()? {
                let id: usize = row.get(0)?;
                let &i = by_id
                    .get(&id)
                    .ok_or_else(|| invalid("event references missing frame"))?;
                let begin: i64 = row.get(1)?;
                let end: i64 = row.get(2)?;
                let mz: f64 = row.get(3)?;
                let width: f64 = row.get(4)?;
                let ce: f64 = row.get(5)?;
                if meta[i].ms_ms_type != expected
                    || begin < 0
                    || end <= begin
                    || end as u64 > meta[i].num_scans as u64
                    || end > u32::MAX as i64
                    || ends.get(&i).is_some_and(|last| begin < *last)
                {
                    return Err(invalid(format!(
                        "invalid or overlapping event at frame {id}"
                    )));
                }
                if !mz.is_finite()
                    || mz <= 0.0
                    || !width.is_finite()
                    || width <= 0.0
                    || !ce.is_finite()
                    || ce < 0.0
                {
                    return Err(invalid(format!("invalid isolation geometry at frame {id}")));
                }
                ends.insert(i, end);
                add(
                    i,
                    begin as u32,
                    end as u32,
                    row.get(6)?,
                    vec![mz.to_bits(), width.to_bits(), ce.to_bits()].into(),
                    msms,
                );
            }
            for (i, f) in meta.iter().enumerate() {
                if f.ms_ms_type == expected && f.num_peaks > 0 && !ends.contains_key(&i) {
                    return Err(invalid(format!(
                        "nonempty frame {} lacks isolation events",
                        f.id
                    )));
                }
            }
        }
        Ok(Some(index))
    }

    pub(crate) fn intervals(&self, frame: usize) -> Vec<(u32, u32)> {
        self.frames[frame]
            .iter()
            .map(|&(g, p)| {
                let o = self.groups[g].observations[p];
                (o.begin, o.end)
            })
            .collect()
    }

    // All run-local read/filter inputs are borrowed; no mutable cache is shared
    // across workers, so streaming and parallel writing use identical evidence.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn keep_mask(
        &self,
        reader: &FrameReader,
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
                    return Err(DnoiseError::Cancelled);
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
                        let f = reader.get(j).map_err(|e| DnoiseError::FrameRead {
                            index: j,
                            message: e.to_string(),
                        })?;
                        e.insert(FlatFrame::from_frame(&f));
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
