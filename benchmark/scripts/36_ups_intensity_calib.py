#!/usr/bin/env python3
"""Calibrate matched per-point intensity thresholds for the UPS intensity-control
arms, at MS1 and at MS1+MS2, to the streak filter's removal on the SAME file.

For each DDA UPS dataset:
  T1 (MS1)  -> ladder of thresholds (with --ms1-polygon, matching the streak arm's
               gate); pick T1 whose kept MS1 fraction matches the streak `denoised`
               arm's kept MS1 fraction.
  T2 (MS2)  -> ladder on the MS2 threshold (--denoise-msms, MS2 collapsed); pick T2
               whose kept MS2 fraction matches the streak `msms` arm's kept MS2.
Then write per-dataset configs:
  dnoise.intensity.<ds>.toml        (MS1-only, T1)         -> arm `intensity`
  dnoise.intensity_msms.<ds>.toml   (MS1 T1 + MS2 T2)      -> arm `intensity_msms`

Usage: uv run scripts/36_ups_intensity_calib.py
"""

from __future__ import annotations

import glob
import shutil
import sqlite3
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DNOISE = (ROOT / ".." / "target" / "release" / "dnoise").resolve()
MS1_BASE = ROOT / "config" / "dnoise.ups_int_ms1.base.toml"
MSMS_BASE = ROOT / "config" / "dnoise.ups_int_msms.base.toml"
DATASETS = ["dda_ups_30spd", "dda_ups_15spd"]
LADDER = [20, 50, 100, 200, 400, 800, 1600, 3200]


def peaks(d: str, ms_level: str) -> int:
    op = "=" if ms_level == "ms1" else "!="
    c = sqlite3.connect(d + "/analysis.tdf")
    v = c.execute(f"SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType{op}0").fetchone()[0]
    c.close()
    return int(v)


def interp(pts, target):
    """pts: list of (T, kept_frac) sorted by T; return T matching target kept_frac."""
    pts = sorted(pts)
    best = min(pts, key=lambda p: abs(p[1] - target))[0]
    for (t0, k0), (t1, k1) in zip(pts, pts[1:]):
        if (k0 - target) * (k1 - target) <= 0 and k0 != k1:
            return round(t0 + (target - k0) * (t1 - t0) / (k1 - k0))
    return best


def _kept(raw, base_cfg, level, flags, T, raw_peaks, tmp):
    out = f"{tmp}/c.d"
    shutil.rmtree(out, ignore_errors=True)
    subprocess.run([str(DNOISE), raw, out, "--config", str(base_cfg), *flags(T)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    k = peaks(out, level) / raw_peaks
    shutil.rmtree(out, ignore_errors=True)
    return k


def calibrate(raw, base_cfg, level, flags, target, raw_peaks, tol=0.01):
    """Coarse ladder to bracket the target kept-fraction, then integer bisection
    (kept decreases monotonically with T) until within `tol` of target. The curve
    is steep, so linear interpolation overshoots -- bisection pins it."""
    with tempfile.TemporaryDirectory() as tmp:
        pts = []
        for T in LADDER:
            k = _kept(raw, base_cfg, level, flags, T, raw_peaks, tmp)
            pts.append((T, k))
            print(f"    T={T:>5d}: keeps {k*100:5.1f}% of {level.upper()}")
        pts.sort()
        # bracket: lo has kept>=target, hi has kept<target
        lo = hi = None
        for (t0, k0), (t1, k1) in zip(pts, pts[1:]):
            if k0 >= target >= k1:
                lo, hi = t0, t1
                break
        if lo is None:  # target outside ladder; fall back to closest
            return min(pts, key=lambda p: abs(p[1] - target))[0]
        best = min(pts, key=lambda p: abs(p[1] - target))
        while hi - lo > 1:
            mid = (lo + hi) // 2
            k = _kept(raw, base_cfg, level, flags, mid, raw_peaks, tmp)
            print(f"    bisect T={mid:>5d}: keeps {k*100:5.1f}% (target {target*100:.1f}%)")
            if abs(k - target) < abs(best[1] - target):
                best = (mid, k)
            if k > target:
                lo = mid
            else:
                hi = mid
            if abs(k - target) <= tol:
                break
        return best[0]


def write_cfg(path, base, repl):
    text = base.read_text()
    import re
    for key, val in repl.items():
        text = re.sub(rf"(?m)^{key}\s*=.*$", f"{key} = {val}", text)
    path.write_text(text)


def main() -> int:
    for ds in DATASETS:
        raw = sorted(glob.glob(f"{ROOT}/data/{ds}/raw/*.d"))[0]
        den = sorted(glob.glob(f"{ROOT}/data/{ds}/denoised/*.d"))[0]
        msms = sorted(glob.glob(f"{ROOT}/data/{ds}/denoised_msms/*.d"))[0]
        r1, r2 = peaks(raw, "ms1"), peaks(raw, "ms2")
        tgt1 = peaks(den, "ms1") / r1                # streak MS1 kept frac
        tgt2 = peaks(msms, "ms2") / r2               # streak MS2 kept frac
        print(f"\n=== {ds}: streak keeps {tgt1*100:.1f}% MS1, {tgt2*100:.1f}% MS2 (targets) ===")
        print("  calibrating T1 (MS1):")
        t1 = calibrate(raw, MS1_BASE, "ms1",
                       lambda T: ["--ms1-polygon", "--min-window-intensity", str(T)], tgt1, r1)
        print("  calibrating T2 (MS2):")
        t2 = calibrate(raw, MSMS_BASE, "ms2",
                       lambda T: ["--ms1-polygon", "--denoise-msms",
                                  "--msms-min-window-intensity", str(T)], tgt2, r2)
        print(f"  -> T1={t1} (MS1), T2={t2} (MS2)")

        write_cfg(ROOT / "config" / f"dnoise.intensity.{ds}.toml", MS1_BASE,
                  {"min_window_intensity": t1})
        write_cfg(ROOT / "config" / f"dnoise.intensity_msms.{ds}.toml", MSMS_BASE,
                  {"min_window_intensity": t1, "msms_min_window_intensity": t2})
        print(f"  wrote dnoise.intensity.{ds}.toml + dnoise.intensity_msms.{ds}.toml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
