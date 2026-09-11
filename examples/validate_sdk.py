"""Independently validate every frame using a user-supplied Bruker SDK library.

Usage: python examples/validate_sdk.py --sdk /path/to/libtimsdata.so output.d
The proprietary SDK is not distributed with dnoise.
"""

import argparse
from contextlib import closing
import ctypes as c
import pathlib
import sqlite3


def validate(lib, path):
    path = path.resolve()
    handle = lib.tims_open(str(path).encode(), 0)
    if not handle:
        raise RuntimeError(sdk_error(lib))
    count = points = 0
    try:
        with closing(sqlite3.connect((path / "analysis.tdf").as_uri() + "?mode=ro", uri=True)) as db:
            for frame, scans, peaks in db.execute(
                "SELECT Id, NumScans, NumPeaks FROM Frames ORDER BY Id"
            ):
                size = 4 * (scans + 2 * peaks)
                buffer = (c.c_uint32 * max(1, size // 4))()
                needed = lib.tims_read_scans_v2(handle, frame, 0, scans, buffer, size)
                if not needed:
                    raise RuntimeError((frame, sdk_error(lib)))
                if needed > size:
                    raise RuntimeError(("buffer too small", frame, needed, size))
                if sum(buffer[:scans]) != peaks:
                    raise RuntimeError(("count mismatch", frame))
                count += 1
                points += peaks
    finally:
        lib.tims_close(handle)
    print(path.name, count, "SDK frames read;", points, "points")


def sdk_error(lib):
    error = c.create_string_buffer(1024)
    lib.tims_get_last_error_string(error, len(error))
    return error.value.decode(errors="replace")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdk", required=True, help="Bruker libtimsdata.so / timsdata.dll")
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args()
    lib = c.CDLL(str(pathlib.Path(args.sdk).resolve()))
    lib.tims_open.argtypes = [c.c_char_p, c.c_uint32]
    lib.tims_open.restype = c.c_uint64
    lib.tims_close.argtypes = [c.c_uint64]
    lib.tims_close.restype = None
    lib.tims_read_scans_v2.argtypes = [
        c.c_uint64, c.c_int64, c.c_uint32, c.c_uint32, c.c_void_p, c.c_uint32
    ]
    lib.tims_read_scans_v2.restype = c.c_uint32
    lib.tims_get_last_error_string.argtypes = [c.c_char_p, c.c_uint32]
    lib.tims_get_last_error_string.restype = c.c_uint32
    for value in args.inputs:
        validate(lib, pathlib.Path(value))


if __name__ == "__main__":
    main()
