"""Search-free verifier regressions: python -m unittest discover -s examples -p 'test_validate_acquisitions.py'."""
import sqlite3
import io
import zstandard as zstd
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
import validate_acquisitions as validation


class ValidationTests(unittest.TestCase):
    def test_subset_and_fragment_identity(self):
        raw = (np.array([1, 2, 4], dtype=np.uint64), np.array([7, 9, 11], dtype=np.uint32))
        subset = (raw[0][[0, 2]], raw[1][[0, 2]])
        validation.compare_points(raw, subset)
        validation.compare_points(raw, raw, identical=True)
        with self.assertRaisesRegex(ValueError, "fragment spectra"):
            validation.compare_points(raw, subset, identical=True)
        for keys, ints in [([1, 5], [7, 11]), ([1, 3], [7, 11]), ([1, 4], [8, 11])]:
            with self.assertRaises(ValueError):
                validation.compare_points(raw, (np.array(keys), np.array(ints)))
        empty = (np.array([], dtype=np.uint64), np.array([], dtype=np.uint32))
        validation.compare_points(empty, empty, identical=True)
        with self.assertRaises(ValueError):
            validation.compare_points(empty, raw)

    def test_metadata_allows_writer_fields_but_rejects_changed_geometry(self):
        with sqlite3.connect(":memory:") as a, sqlite3.connect(":memory:") as b:
            for db in (a, b):
                db.executescript("CREATE TABLE Frames (Id INTEGER, TimsId INTEGER, NumPeaks INTEGER, Time REAL);"
                                 "INSERT INTO Frames VALUES (1,64,2,1.0);"
                                 "CREATE TABLE Windows (Mz REAL); INSERT INTO Windows VALUES (500),(600);")
            b.execute("UPDATE Frames SET TimsId=100, NumPeaks=1")
            validation.compare_metadata(a, b)
            b.executescript("DELETE FROM Windows; INSERT INTO Windows VALUES (600),(500)")
            validation.compare_metadata(a, b)
            b.execute("UPDATE Windows SET Mz=700 WHERE Mz=600")
            with self.assertRaisesRegex(ValueError, "Windows"):
                validation.compare_metadata(a, b)
            b.execute("UPDATE Windows SET Mz=600 WHERE Mz=700")
            b.execute("UPDATE Frames SET Time=2.0")
            with self.assertRaisesRegex(ValueError, "timing"):
                validation.compare_metadata(a, b)

    def test_report_is_checked_against_sdk_totals(self):
        raw = {"ms1_intensity": 10, "msms_intensity": 20, "ms1_points": 2, "msms_points": 3}
        output = {"ms1_intensity": 5, "msms_intensity": 0, "ms1_points": 1, "msms_points": 0}
        stats = {"raw_ms1_summed_intensity": 10, "raw_msms_summed_intensity": 20,
                 "kept_ms1_summed_intensity": 5, "kept_msms_summed_intensity": 0,
                 "raw_points": 5, "kept_points": 1}
        validation.check_report(stats, raw, output)
        stats["raw_ms1_summed_intensity"] = 9
        with self.assertRaisesRegex(ValueError, "disagrees"):
            validation.check_report(stats, raw, output)

    def test_native_decoder_preserves_stored_integers_and_empty_scans(self):
        values = np.array([3, 4, 0, 101, 7718, 4, 20, 201, 30], dtype="<u4")
        payload = zstd.ZstdCompressor().compress(values.view(np.uint8).reshape(-1, 4).T.tobytes())
        record = (len(payload) + 8).to_bytes(4, "little") + (3).to_bytes(4, "little") + payload
        keys, intensity = validation.read_native(io.BytesIO(record), 0, (1, 3, 3, 0))
        np.testing.assert_array_equal(keys, [100, 104, (2 << 32) + 200])
        np.testing.assert_array_equal(intensity, [7718, 20, 30])
        with self.assertRaisesRegex(ValueError, "truncated"):
            validation.read_native(io.BytesIO(record[:-1]), 0, (1, 3, 3, 0))
        with self.assertRaisesRegex(ValueError, "size mismatch"):
            validation.read_native(io.BytesIO(record), 0, (1, 3, 4, 0))
        empty = (8).to_bytes(4, "little") + (3).to_bytes(4, "little")
        self.assertEqual(len(validation.read_native(io.BytesIO(empty), 0, (1, 3, 0, 0))[0]), 0)

    def test_existing_output_is_rejected_before_sdk_or_processing(self):
        with tempfile.TemporaryDirectory() as folder:
            sentinel = Path(folder) / "sentinel"
            sentinel.write_text("keep")
            with patch("sys.argv", ["validate", "--binary", "unused", "--sdk", "unused",
                                    "--output-root", folder]):
                with self.assertRaisesRegex(ValueError, "already exist"):
                    validation.main()
            self.assertEqual(sentinel.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
