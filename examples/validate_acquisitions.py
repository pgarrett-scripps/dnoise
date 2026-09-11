"""Run MS1, MS/MS, and neighbor recipes and independently check every frame.

Requires Python 3.10+, numpy, zstandard, a built dnoise binary, and a user-supplied Bruker
SDK. Downloads nothing. Existing output roots are never reused or overwritten.
See docs/acquisition-validation.md for the manifest and result definitions.
"""
import argparse
import ctypes as c
from contextlib import ExitStack
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3
import subprocess
import sys

import numpy as np
import zstandard as zstd

NATIVE_FILES = ("analysis.tdf", "analysis.tdf_bin")
WRITER_COLUMNS = {"TimsId", "NumPeaks", "MaxIntensity", "SummedIntensities"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fingerprints(path):
    return {name: {"bytes": (path / name).stat().st_size,
                   "sha256": sha256(path / name)} for name in NATIVE_FILES}


def connect(path):
    return sqlite3.connect((path / "analysis.tdf").resolve().as_uri() + "?mode=ro", uri=True)


def same_rows(a, b, query, description):
    sentinel = object()
    require(all(x == y for x, y in itertools.zip_longest(
        a.execute(query), b.execute(query), fillvalue=sentinel)), description)


def compare_metadata(a, b):
    same_rows(a, b, "SELECT type,name,sql FROM sqlite_master ORDER BY type,name",
              "SQLite schema changed")
    tables = [r[0] for r in a.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name!='Frames'")]
    for table in tables:
        quoted = '"' + table.replace('"', '""') + '"'
        # Sort full rows so equivalent tables do not depend on physical row order.
        count = len(a.execute(f"SELECT * FROM {quoted} LIMIT 0").description)
        order = ','.join(str(i + 1) for i in range(count))
        same_rows(a, b, f"SELECT * FROM {quoted} ORDER BY {order}",
                  f"Metadata table changed: {table}")
    columns = [r[1] for r in a.execute("PRAGMA table_info(Frames)")
               if r[1] not in WRITER_COLUMNS]
    fields = ','.join('"' + x.replace('"', '""') + '"' for x in columns)
    same_rows(a, b, f"SELECT {fields} FROM Frames ORDER BY Id",
              "Frame identities, timing or acquisition metadata changed")
    return len(tables)


class SDK:
    def __init__(self, path):
        self.lib = c.CDLL(str(path))
        self.lib.tims_open.argtypes = [c.c_char_p, c.c_uint32]
        self.lib.tims_open.restype = c.c_uint64
        self.lib.tims_close.argtypes = [c.c_uint64]
        self.lib.tims_close.restype = None
        self.lib.tims_read_scans_v2.argtypes = [
            c.c_uint64, c.c_int64, c.c_uint32, c.c_uint32, c.c_void_p, c.c_uint32]
        self.lib.tims_read_scans_v2.restype = c.c_uint32
        self.lib.tims_get_last_error_string.argtypes = [c.c_char_p, c.c_uint32]
        self.lib.tims_get_last_error_string.restype = c.c_uint32

    def error(self):
        buffer = c.create_string_buffer(2048)
        self.lib.tims_get_last_error_string(buffer, len(buffer))
        return buffer.value.decode(errors="replace")

    def open(self, path, stack):
        handle = self.lib.tims_open(str(path).encode(), 0)
        require(handle != 0, f"SDK could not open {path}: {self.error()}")
        stack.callback(self.lib.tims_close, handle)
        return handle

    def decode(self, handle, row):
        frame, scans, peaks, _ = row
        size = 4 * (scans + 2 * peaks)
        require(0 <= scans <= 1000000 and 0 <= peaks and size <= 512 * 1024 * 1024,
                f"Frame {frame}: invalid or excessive decoded allocation")
        if scans == 0:
            require(peaks == 0, f"Frame {frame}: peaks without scans")
            return np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.uint32)
        buffer = (c.c_uint32 * (size // 4))()
        needed = self.lib.tims_read_scans_v2(handle, frame, 0, scans, buffer, size)
        require(0 < needed <= size, f"Frame {frame}: SDK read failed: {self.error()}")
        values = np.frombuffer(buffer, dtype=np.uint32)
        counts = values[:scans]
        require(int(counts.sum(dtype=np.uint64)) == peaks, f"Frame {frame}: SDK count mismatch")
        keys = np.empty(peaks, dtype=np.uint64)
        intensities = np.empty(peaks, dtype=np.uint32)
        pos, out = scans, 0
        for scan, count in enumerate(counts):
            n = int(count)
            if n:
                keys[out:out+n] = values[pos:pos+n].astype(np.uint64) + (scan << 32)
                intensities[out:out+n] = values[pos+n:pos+2*n]
            pos += 2 * n
            out += n
        require(bool(np.all(keys[1:] > keys[:-1])), f"Frame {frame}: unordered or duplicate points")
        return keys, intensities


def read_native(file, offset, row):
    """Independent, bounded type-2 decoding for exact stored-integer QC totals."""
    frame, scans, peaks, _ = row
    file.seek(offset)
    header = file.read(8)
    require(len(header) == 8, f"Frame {frame}: short native header")
    size, header_scans = int.from_bytes(header[:4], "little"), int.from_bytes(header[4:], "little")
    require(header_scans == scans and 8 <= size <= 256 * 1024 * 1024,
            f"Frame {frame}: invalid native header")
    if size == 8:
        require(peaks == 0, f"Frame {frame}: missing native peaks")
        return np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.uint32)
    payload = file.read(size - 8)
    require(len(payload) == size - 8, f"Frame {frame}: truncated native record")
    expected = 4 * (scans + 2 * peaks)
    require(0 < expected <= 256 * 1024 * 1024, "Excessive native decoded size")
    require(zstd.frame_content_size(payload) == expected, "Native payload size mismatch")
    decoded = zstd.ZstdDecompressor().decompress(payload, max_output_size=expected)
    values = np.frombuffer(decoded, dtype=np.uint8).reshape(4, -1).T.copy().view("<u4").reshape(-1)
    require(scans > 0 and int(values[0]) == scans, "Native scan header mismatch")
    counts = np.empty(scans, dtype=np.int64)
    require(bool(np.all(values[1:scans] % 2 == 0)), "Invalid native scan counts")
    counts[:-1] = values[1:scans] // 2
    counts[-1] = peaks - int(counts[:-1].sum())
    require(bool(np.all(counts >= 0)), "Native scan counts exceed point count")
    deltas, intensities = values[scans::2], values[scans+1::2]
    keys = np.empty(peaks, dtype=np.uint64)
    pos = 0
    for scan, count in enumerate(counts):
        n = int(count)
        if n:
            offsets = deltas[pos:pos+n].cumsum(dtype=np.uint64)
            require(bool(np.all(deltas[pos:pos+n] > 0)) and int(offsets[-1]) <= 2**32,
                    "Invalid native TOF deltas")
            keys[pos:pos+n] = (scan << 32) + offsets - 1
        pos += n
    return keys, intensities


def compare_points(raw, output, identical=False):
    raw_keys, raw_intensity = raw
    keys, intensity = output
    positions = np.searchsorted(raw_keys, keys)
    require(bool(np.all(positions < len(raw_keys))), "Output added native coordinates")
    require(np.array_equal(raw_keys[positions], keys), "Output changed native coordinates")
    require(np.array_equal(raw_intensity[positions], intensity), "Output changed native intensities")
    if identical:
        require(np.array_equal(raw_keys, keys), "MS1-only recipe changed fragment spectra")


def verify(sdk, source, outputs):
    paths = [source] + [p for _, p in outputs]
    names = ["raw"] + [name for name, _ in outputs]
    totals = {name: {"ms1_points": 0, "msms_points": 0,
                     "ms1_intensity": 0, "msms_intensity": 0} for name in names}
    sdk_totals = {name: dict(value) for name, value in totals.items()}
    with ExitStack() as stack:
        dbs = [connect(p) for p in paths]
        for db in dbs:
            stack.callback(db.close)
        handles = [sdk.open(p, stack) for p in paths]
        native_files = [stack.enter_context((p / "analysis.tdf_bin").open("rb")) for p in paths]
        offsets = [dict(db.execute("SELECT Id,TimsId FROM Frames")) for db in dbs]
        rows = [list(db.execute("SELECT Id,NumScans,NumPeaks,MsMsType FROM Frames ORDER BY Id"))
                for db in dbs]
        for j in range(1, len(paths)):
            require(len(rows[0]) == len(rows[j]), "Frame count changed")
            compare_metadata(dbs[0], dbs[j])
        for i, row in enumerate(rows[0]):
            raw = sdk.decode(handles[0], row)
            raw_native = read_native(native_files[0], offsets[0][row[0]], row)
            for j, name in enumerate(names):
                f, scans, peaks, level = rows[j][i]
                require((f, scans, level) == (row[0], row[1], row[3]), "Frame metadata mismatch")
                points = raw if j == 0 else sdk.decode(handles[j], rows[j][i])
                if j:
                    compare_points(raw, points, identical=name == "ms1" and level != 0)
                native = raw_native if j == 0 else read_native(native_files[j], offsets[j][f], rows[j][i])
                require(np.array_equal(points[0], native[0]), "SDK and native coordinates disagree")
                if j:
                    compare_points(raw_native, native, identical=name == "ms1" and level != 0)
                prefix = "ms1" if level == 0 else "msms"
                totals[name][prefix + "_points"] += peaks
                totals[name][prefix + "_intensity"] += int(native[1].sum(dtype=np.uint64))
                sdk_totals[name][prefix + "_points"] += peaks
                sdk_totals[name][prefix + "_intensity"] += int(points[1].sum(dtype=np.uint64))
            if i % 1000 == 0:
                print(f"{source.name}: SDK verified {i}/{len(rows[0])} frames", flush=True)
    return {"frames_per_file": len(rows[0]), "default_msms_identical": True,
            "all_outputs_native_subsets": True, "metadata_preserved": True, "intensity_basis": "stored_integer",
            "totals": totals, "sdk_totals": sdk_totals}


def check_report(stats, raw, output):
    for level in ("ms1", "msms"):
        require(stats[f"raw_{level}_summed_intensity"] == raw[level + "_intensity"],
                f"{level} input intensity report disagrees with independent native decoding")
        require(stats[f"kept_{level}_summed_intensity"] == output[level + "_intensity"],
                f"{level} output intensity report disagrees with independent native decoding")
    require(stats["raw_points"] == raw["ms1_points"] + raw["msms_points"], "Input point report mismatch")
    require(stats["kept_points"] == output["ms1_points"] + output["msms_points"], "Output point report mismatch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path(__file__).resolve().parents[1] / "test-data/acquisition-examples.json")
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--sdk", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True, help="New directory for all outputs and evidence")
    parser.add_argument("--method", action="append", help="Select exact manifest method(s); default all")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--frame-batch-size", type=int, default=32)
    args = parser.parse_args()
    require(args.threads > 0 and args.frame_batch_size > 0, "Threads and batch size must be positive")
    manifest_path, binary, sdk_path = (x.resolve() for x in (args.manifest, args.binary, args.sdk))
    root = args.output_root.resolve()
    require(not root.exists(), "Output root must not already exist")
    manifest = json.loads(manifest_path.read_text())
    require(manifest["schema_version"] == 1, "Unsupported manifest version")
    examples = [e for e in manifest["examples"] if not args.method or e["method"] in args.method]
    require(bool(examples), "No examples selected")
    require(not args.method or set(args.method) <= {e["method"] for e in examples}, "Unknown method selected")
    sdk = SDK(sdk_path)
    selected = []
    for entry in examples:
        source = (manifest_path.parent / entry["directory"]).resolve()
        require(not source.is_relative_to(root) and not root.is_relative_to(source), "Input and output paths overlap")
        actual = fingerprints(source)
        require(actual == {name: entry["files"][name] for name in NATIVE_FILES}, f"Source checksum mismatch: {source}")
        with ExitStack() as stack:
            db = connect(source)
            stack.callback(db.close)
            require(db.execute("SELECT Value FROM GlobalMetadata WHERE Key='TimsCompressionType'").fetchone()[0] == "2", "Only compression type 2 is supported")
            kinds = {r[0] for r in db.execute("SELECT DISTINCT MsMsType FROM Frames") if r[0] != 0}
        require(len(kinds) <= 1 and kinds <= {8, 9, 10}, f"Unsupported/mixed acquisition: {kinds}")
        selected.append((entry, source, actual, next(iter(kinds), 0)))
    binary_hash, sdk_hash = sha256(binary), sha256(sdk_path)
    root.mkdir(parents=True)
    result = {"schema_version": 1, "status": "running", "binary": str(binary),
              "binary_sha256": binary_hash, "sdk_sha256": sdk_hash,
              "manifest_sha256": sha256(manifest_path), "examples": [],
              "verifier_sha256": sha256(Path(__file__)), "python_version": sys.version,
              "numpy_version": np.__version__, "zstandard_version": zstd.__version__}
    try:
        for index, (entry, source, original, kind) in enumerate(selected):
            item = {"method": entry["method"], "input": str(source), "input_files": original,
                    "recipes": {}}
            result["examples"].append(item)
            folder = root / f"example-{index+1}"
            folder.mkdir()
            outputs = []
            for name in ("ms1", "msms", "neighbors"):
                output = folder / (name + ".d")
                report = folder / (name + ".json")
                command = [str(binary), str(source), str(output), "--threads", str(args.threads),
                           "--frame-batch-size", str(args.frame_batch_size), "--report", str(report)]
                if name != "ms1":
                    command += ["--denoise-msms"]
                if name == "neighbors":
                    command += ["--ms1-neighbor-radius", "1"]
                    if kind in (9, 10):
                        command += ["--dia-neighbor-radius" if kind == 9 else "--prm-neighbor-radius", "1"]
                item["recipes"][name] = {"command": command}
                print(f"{entry['method']}: processing {name}", flush=True)
                with (folder / (name + ".log")).open("w") as log:
                    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
                item["recipes"][name]["report"] = json.loads(report.read_text())
                item["recipes"][name]["output_files"] = fingerprints(output)
                outputs.append((name, output))
            item["verification"] = verify(sdk, source, outputs)
            for name, _ in outputs:
                check_report(item["recipes"][name]["report"]["stats"],
                             item["verification"]["totals"]["raw"],
                             item["verification"]["totals"][name])
            require(fingerprints(source) == original, "Source files changed during validation")
        require(sha256(binary) == binary_hash and sha256(sdk_path) == sdk_hash,
                "Binary or SDK changed during validation")
        result["status"] = "passed"
    except BaseException as error:
        result["status"] = "failed"
        result["error"] = str(error)
        raise
    finally:
        (root / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"Passed: {root / 'validation.json'}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
