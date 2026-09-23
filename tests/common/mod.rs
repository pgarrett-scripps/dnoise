//! Synthetic fixtures generated from scratch; no vendor data redistributed.
#![allow(dead_code)]
use rusqlite::{Connection, params};
use std::io::Write;
use std::path::{Path, PathBuf};
pub fn fixture(root: &Path, name: &str, dia: bool) -> PathBuf {
    let path = root.join(name);
    std::fs::create_dir(&path).unwrap();
    let db = Connection::open(path.join("analysis.tdf")).unwrap();
    db.execute_batch("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT);
        CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL, ScanMode INTEGER, MsMsType INTEGER, TimsId INTEGER,
          MaxIntensity INTEGER, SummedIntensities INTEGER, NumScans INTEGER, NumPeaks INTEGER, AccumulationTime REAL,
          MzCalibration INTEGER, TimsCalibration INTEGER);
        CREATE TABLE TimsCalibration (Id INTEGER PRIMARY KEY, ModelType INTEGER, C0 REAL, C1 REAL,
          C2 REAL, C3 REAL, C4 REAL, C5 REAL, C6 REAL, C7 REAL, C8 REAL, C9 REAL);
        INSERT INTO TimsCalibration VALUES (1, 2, 1, 15, 239.3, 103.0, 33.6, 1, -0.0266, 171.4, 16.8, 1732.7);").unwrap();
    for (key, value) in [
        ("TimsCompressionType", "2"),
        ("AcquisitionSoftware", "synthetic fixture"),
        ("MzAcqRangeLower", "100"),
        ("MzAcqRangeUpper", "1700"),
        ("DigitizerNumSamples", "100000"),
        ("OneOverK0AcqRangeLower", "0.5"),
        ("OneOverK0AcqRangeUpper", "1.5"),
    ] {
        db.execute(
            "INSERT INTO GlobalMetadata VALUES (?1,?2)",
            params![key, value],
        )
        .unwrap();
    }
    if dia {
        db.execute_batch("CREATE TABLE DiaFrameMsMsInfo (Frame INTEGER, WindowGroup INTEGER);
          CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, ScanNumEnd INTEGER,
            IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL);
          INSERT INTO DiaFrameMsMsInfo VALUES (2,1),(4,1);
          INSERT INTO DiaFrameMsMsWindows VALUES (1,0,16,900,1600,30);").unwrap();
    }
    let mut bin = std::fs::File::create(path.join("analysis.tdf_bin")).unwrap();
    bin.write_all(&[0u8; 64]).unwrap();
    let mut offset = 64;
    for id in 1..=4 {
        let mut points: Vec<_> = (3..11).map(|s| (s, 10000, 100)).collect();
        points.push((1, 60000, 1));
        if id == 4 {
            points.clear();
        }
        let record = if points.is_empty() {
            dnoise::codec::encode_empty_frame_type2(16)
        } else {
            dnoise::codec::encode_frame_type2(16, &points)
        };
        db.execute(
            "INSERT INTO Frames VALUES (?1,?2,8,?3,?4,100,?5,16,?6,100,1,1)",
            params![
                id,
                id as f64 * 10.,
                if id % 2 == 0 {
                    if dia { 9 } else { 8 }
                } else {
                    0
                },
                offset,
                points.iter().map(|p| p.2 as u64).sum::<u64>(),
                points.len()
            ],
        )
        .unwrap();
        bin.write_all(&record).unwrap();
        offset += record.len();
    }
    path
}

/// PRM schema matching PXD049405, including touching events and an empty frame.
pub fn prm_fixture(root: &Path, name: &str) -> PathBuf {
    let path = fixture(root, name, false);
    let db = Connection::open(path.join("analysis.tdf")).unwrap();
    db.execute_batch("UPDATE Frames SET ScanMode=10, MsMsType=10 WHERE MsMsType=8;
        CREATE TABLE PrmTargets (Id INTEGER PRIMARY KEY, ExternalId TEXT, Time REAL NOT NULL,
            OneOverK0 REAL NOT NULL, MonoisotopicMz REAL NOT NULL, Charge INTEGER NOT NULL, Description TEXT NOT NULL);
        INSERT INTO PrmTargets VALUES (1,NULL,20,0.8,500,2,'target A'), (1000001,'heavy',20,1.0,600,2,'target B');
        CREATE TABLE PrmFrameMsMsInfo (Frame INTEGER NOT NULL, ScanNumBegin INTEGER NOT NULL,
            ScanNumEnd INTEGER NOT NULL, IsolationMz REAL NOT NULL, IsolationWidth REAL NOT NULL,
            CollisionEnergy REAL NOT NULL, Target INTEGER NOT NULL, PRIMARY KEY(Frame, ScanNumBegin),
            FOREIGN KEY(Frame) REFERENCES Frames(Id), FOREIGN KEY(Target) REFERENCES PrmTargets(Id)) WITHOUT ROWID;
        INSERT INTO PrmFrameMsMsInfo VALUES (2,3,7,500,2,30,1),(2,7,11,600,2,30,1000001),(4,3,11,500,2,30,1);
        CREATE TABLE PrmFrameMeasurementMode (Frame INTEGER PRIMARY KEY, MeasurementModeId TEXT);
        INSERT INTO PrmFrameMeasurementMode VALUES (2,NULL),(4,NULL);
        CREATE TABLE DiaFrameMsMsInfo (Frame INTEGER PRIMARY KEY, WindowGroup INTEGER);
        CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, ScanNumEnd INTEGER,
            IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL);").unwrap();
    path
}
