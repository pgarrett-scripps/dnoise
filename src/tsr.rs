//! Thin adapter over `timsrust` 0.6.
//!
//! timsrust 0.6 split into sub-crates and changed the reader API: frames are
//! indexed by `Frames.Id` (1-based) instead of position, ion arrays are typed
//! newtypes behind accessors, and the metadata no longer carries the m/z and
//! 1/K0 converters. Its converters also take integer indices and truncate on
//! the inverse, and the linear 1/K0 line ends at `max(NumScans) - 1` instead
//! of `max(NumScans)`. dnoise needs fractional inverses and must reproduce
//! earlier output, so this module keeps the 0.4 shapes and formulas:
//!
//! - [`FrameReader::get`] takes the 0-based position (`Frames.Id - 1` order);
//! - [`Frame`] exposes plain `scan_offsets` / `tof_indices` / `intensities`;
//! - [`Tof2MzConverter`] and [`Scan2ImConverter`] are timsrust 0.4.2's
//!   straight lines, built from the same `GlobalMetadata` keys.

use std::collections::HashMap;
use std::path::Path;

use rusqlite::{Connection, OpenFlags};

/// Scalar conversion with a fractional inverse (timsrust 0.4's trait).
pub trait ConvertableDomain {
    /// Forward conversion (index -> physical value).
    fn convert<T: Into<f64> + Copy>(&self, value: T) -> f64;
    /// Inverse conversion (physical value -> fractional index).
    fn invert<T: Into<f64> + Copy>(&self, value: T) -> f64;
}

/// TOF index -> m/z: `mz = (sqrt(mz_min) + slope * tof)^2` (timsrust 0.4.2).
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Tof2MzConverter {
    tof_intercept: f64,
    tof_slope: f64,
}

impl Tof2MzConverter {
    /// Line between `mz_min` at TOF 0 and `mz_max` at `tof_max_index`.
    pub fn from_boundaries(mz_min: f64, mz_max: f64, tof_max_index: u32) -> Self {
        let tof_intercept = mz_min.sqrt();
        let tof_slope = (mz_max.sqrt() - tof_intercept) / tof_max_index as f64;
        Self {
            tof_intercept,
            tof_slope,
        }
    }
}

impl ConvertableDomain for Tof2MzConverter {
    fn convert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        let tof_index: f64 = value.into();
        (self.tof_intercept + self.tof_slope * tof_index).powi(2)
    }
    fn invert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        let mz_value: f64 = value.into();
        (mz_value.sqrt() - self.tof_intercept) / self.tof_slope
    }
}

/// Scan index -> 1/K0: straight line from `im_max` at scan 0 to `im_min` at
/// `max(NumScans)` (timsrust 0.4.2).
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub struct Scan2ImConverter {
    scan_intercept: f64,
    scan_slope: f64,
}

impl Scan2ImConverter {
    /// Line between `im_max` at scan 0 and `im_min` at `scan_max_index`.
    pub fn from_boundaries(im_min: f64, im_max: f64, scan_max_index: u32) -> Self {
        let scan_intercept = im_max;
        let scan_slope = (im_min - scan_intercept) / scan_max_index as f64;
        Self {
            scan_intercept,
            scan_slope,
        }
    }
}

impl ConvertableDomain for Scan2ImConverter {
    fn convert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        let scan_index: f64 = value.into();
        self.scan_intercept + self.scan_slope * scan_index
    }
    fn invert<T: Into<f64> + Copy>(&self, value: T) -> f64 {
        let im_value: f64 = value.into();
        (im_value - self.scan_intercept) / self.scan_slope
    }
}

/// The run-level converters dnoise reads from `analysis.tdf`.
#[derive(Debug, Clone, Copy)]
pub struct Metadata {
    /// TOF index -> m/z.
    pub mz_converter: Tof2MzConverter,
    /// Scan index -> 1/K0 (uncalibrated linear scale).
    pub im_converter: Scan2ImConverter,
}

/// Error reading the run metadata.
#[derive(Debug, thiserror::Error)]
pub enum MetadataReaderError {
    /// SQLite error.
    #[error("{0}")]
    Sql(#[from] rusqlite::Error),
    /// A required `GlobalMetadata` key is missing.
    #[error("Key not found: {0}")]
    KeyNotFound(String),
    /// A `GlobalMetadata` value does not parse.
    #[error("Key not parsable: {0}")]
    ParseError(String),
    /// The `Frames` table is empty.
    #[error("Frames table is empty")]
    NoFrames,
}

/// Reads [`Metadata`] the way timsrust 0.4.2's `MetadataReader` did.
pub struct MetadataReader;

const OTOF_CONTROL: &str = "Bruker otofControl";

impl MetadataReader {
    /// `path` is the `.d` folder or its `analysis.tdf`.
    #[allow(clippy::new_ret_no_self)]
    pub fn new(path: impl AsRef<Path>) -> Result<Metadata, MetadataReaderError> {
        let path = path.as_ref();
        let tdf = if path.is_dir() {
            path.join("analysis.tdf")
        } else {
            path.to_path_buf()
        };
        let conn = Connection::open_with_flags(&tdf, OpenFlags::SQLITE_OPEN_READ_ONLY)?;
        let mut stmt = conn.prepare("SELECT Key, Value FROM GlobalMetadata")?;
        // Values are read as TEXT, as 0.4.2 did (a non-text value is an error).
        let kv: HashMap<String, String> = stmt
            .query_map([], |r| Ok((r.get::<_, String>(0)?, r.get::<_, String>(1)?)))?
            .collect::<Result<_, _>>()?;
        fn parse<T: std::str::FromStr>(
            kv: &HashMap<String, String>,
            key: &str,
        ) -> Result<T, MetadataReaderError> {
            kv.get(key)
                .ok_or_else(|| MetadataReaderError::KeyNotFound(key.to_string()))?
                .parse()
                .map_err(|_| MetadataReaderError::ParseError(key.to_string()))
        }
        // 0.4.2 parsed the compression type here and failed when it was missing.
        let _compression_type: u8 = parse(&kv, "TimsCompressionType")?;
        let software = kv
            .get("AcquisitionSoftware")
            .ok_or_else(|| MetadataReaderError::KeyNotFound("AcquisitionSoftware".into()))?;
        let mut mz_min: f64 = parse(&kv, "MzAcqRangeLower")?;
        let mut mz_max: f64 = parse(&kv, "MzAcqRangeUpper")?;
        if software == OTOF_CONTROL {
            mz_min -= 5.0;
            mz_max += 5.0;
        }
        let tof_max_index: u32 = parse(&kv, "DigitizerNumSamples")?;
        let im_min: f64 = parse(&kv, "OneOverK0AcqRangeLower")?;
        let im_max: f64 = parse(&kv, "OneOverK0AcqRangeUpper")?;
        let scan_max: Option<i64> =
            conn.query_row("SELECT MAX(NumScans) FROM Frames", [], |r| r.get(0))?;
        let scan_max = scan_max.ok_or(MetadataReaderError::NoFrames)?;
        Ok(Metadata {
            mz_converter: Tof2MzConverter::from_boundaries(mz_min, mz_max, tof_max_index),
            im_converter: Scan2ImConverter::from_boundaries(im_min, im_max, scan_max as u32),
        })
    }
}

/// One decoded frame in timsrust 0.4's CSR layout.
#[derive(Debug, Clone, Default)]
pub struct Frame {
    /// `Frames.Id`.
    pub index: usize,
    /// CSR row pointer: scan `s` holds points `scan_offsets[s]..scan_offsets[s + 1]`.
    pub scan_offsets: Vec<usize>,
    /// Per-point TOF index.
    pub tof_indices: Vec<u32>,
    /// Per-point raw intensity.
    pub intensities: Vec<u32>,
}

/// Error opening or decoding frames.
#[derive(Debug, thiserror::Error)]
#[error("{0}")]
pub struct FrameReaderError(String);

/// Frame reader indexed by 0-based position in `Frames.Id` order.
pub struct FrameReader {
    inner: timsrust_tdf::TdfFrameReader,
    ids: Vec<usize>,
}

impl std::fmt::Debug for FrameReader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("FrameReader")
            .field("frames", &self.ids.len())
            .finish()
    }
}

impl FrameReader {
    /// Open the `.d` folder `path`.
    ///
    /// Like timsrust 0.4.2, only `TimsCompressionType` 2 is accepted. The
    /// compression type is read here directly so timsrust's `Metadata::new`
    /// (which builds per-frame info) is not needed.
    pub fn new(path: impl AsRef<Path>) -> Result<Self, FrameReaderError> {
        let p = path.as_ref();
        let s = p
            .to_str()
            .ok_or_else(|| FrameReaderError(format!("non UTF-8 path {}", p.display())))?;
        let tdf = if p.is_dir() {
            p.join("analysis.tdf")
        } else {
            p.to_path_buf()
        };
        let compression_type = read_compression_type(&tdf)?;
        if compression_type != 2 {
            return Err(FrameReaderError(format!(
                "Compression type {compression_type} not understood"
            )));
        }
        let inner = timsrust_tdf::TdfFrameReader::without_metadata(s, compression_type, 0)
            .map_err(|e| FrameReaderError(e.to_string()))?;
        // `iter_indices` yields `Frames.Id` in HashMap order; sort so position
        // `i` is the i-th frame by Id, which is 0.4.2's `get(i)` order.
        let mut ids: Vec<usize> = inner.iter_indices().collect();
        ids.sort_unstable();
        Ok(Self { inner, ids })
    }

    /// Number of frames.
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// True when the run has no frames.
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// Decode the frame at 0-based position `i`.
    pub fn get(&self, i: usize) -> Result<Frame, FrameReaderError> {
        let id = *self
            .ids
            .get(i)
            .ok_or_else(|| FrameReaderError(format!("frame index {i} out of bounds")))?;
        let ions = self
            .inner
            .get_ions(id)
            .map_err(|e| FrameReaderError(e.to_string()))?;
        Ok(Frame {
            index: id,
            scan_offsets: ions.scan_offsets().clone(),
            tof_indices: ions.tof_indices().iter().map(|&t| u32::from(t)).collect(),
            intensities: ions.intensities().iter().map(|&x| u32::from(x)).collect(),
        })
    }
}

/// `GlobalMetadata.TimsCompressionType`, parsed as 0.4.2 did (TEXT -> `u8`).
fn read_compression_type(tdf: &Path) -> Result<u8, FrameReaderError> {
    let err = |e: &dyn std::fmt::Display| FrameReaderError(e.to_string());
    let conn =
        Connection::open_with_flags(tdf, OpenFlags::SQLITE_OPEN_READ_ONLY).map_err(|e| err(&e))?;
    let value: String = conn
        .query_row(
            "SELECT Value FROM GlobalMetadata WHERE Key='TimsCompressionType'",
            [],
            |r| r.get(0),
        )
        .map_err(|_| err(&"Key not found: TimsCompressionType"))?;
    value
        .parse()
        .map_err(|_| err(&"Key not parsable: TimsCompressionType"))
}
