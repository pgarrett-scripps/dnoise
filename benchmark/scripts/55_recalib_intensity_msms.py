#!/usr/bin/env python3
"""Calibrate a matched per-point MS/MS intensity threshold (T2) for the main HYE
ddaPASEF benchmark, extending the existing MS1-only intensity control
(dnoise.intensity.toml / .15min.toml) to the MS1+MS/MS level.

T1 (MS1) is NOT recalibrated here -- it is reused as-is from the existing
dnoise.intensity.toml / .15min.toml (86 @ 5min, 82 @ 15min; see
35_recalib_intensity.py). T1 and T2 are independent: MS1 and MS/MS frames are
disjoint (Frames.MsMsType = 0 vs != 0), so the MS1 threshold has no effect on
the MS2 point count used to calibrate T2, and vice versa (same independence
36_ups_intensity_calib.py relies on for the UPS arms).

For each DDA gradient: rebuild a temporary streak MS1+MS/MS copy of REP1 (the
same config/flags 04_denoise.sh uses for the `msms` arm) to measure the streak
filter's MS2-kept-fraction, delete it, then bisect a per-point MS2 threshold T2
(via degenerate CLI flags: msms_mz_half_width=0, msms_min_feature_length=1,
msms_max_internal_gap=0, msms_iterations=1) until its MS2-kept-fraction matches.
Writes T1 (reused) and T2 (new) into config/dnoise.intensity_msms.toml (5min)
and dnoise.intensity_msms.15min.toml (15min).

Storage note: this script only ever materializes ONE temporary single-file .d
copy at a time (in a tempdir), deleted immediately after each measurement --
never the full 18-run arm.

Usage: uv run scripts/55_recalib_intensity_msms.py
"""

from __future__ import annotations

import glob
import re
import shutil
import sqlite3
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DNOISE = (ROOT / ".." / "target" / "release" / "dnoise").resolve()
STREAK_CFG = ROOT / "config" / "dnoise.toml"

# gradient -> (existing MS1-only intensity config (source of T1), output config)
GRADS = {
    "dda_5min": {
        "ms1_cfg": ROOT / "config" / "dnoise.intensity.toml",
        "out_cfg": ROOT / "config" / "dnoise.intensity_msms.toml",
    },
    "dda_15min": {
        "ms1_cfg": ROOT / "config" / "dnoise.intensity.15min.toml",
        "out_cfg": ROOT / "config" / "dnoise.intensity_msms.15min.toml",
    },
}
LADDER = [20, 50, 100, 200, 400, 800, 1600, 3200]

TEMPLATE = """\
# Strict intensity-threshold baseline, MS1+MS/MS (control arm for the streak
# filter's `msms` arm). Both frame types are collapsed to a pure per-point
# threshold: a point is kept iff its own intensity >= the relevant floor. It is
# NOT the streak filter. See dnoise.intensity{suffix}.toml for the MS1-only
# rationale (mz_half_width=0 etc. degenerating the vertical filter to a
# per-point test); the same trick is applied to the MS/MS filter here via
# msms_mz_half_width=0, msms_min_feature_length=1, msms_max_internal_gap=0,
# msms_iterations=1.
#
# T1 = {t1} is reused unchanged from dnoise.intensity{suffix}.toml (calibrated
# by 35_recalib_intensity.py against the streak MS1 arm's removal on this
# gradient). T2 = {t2} was calibrated by 55_recalib_intensity_msms.py to match
# the streak `msms` arm's MS2-point removal on this gradient.

mz_half_width = 0
min_feature_length = 1
max_internal_gap = 0
min_window_intensity = {t1}
min_feature_intensity = 0
iterations = 1
halo = false
denoise_msms = true
msms_mz_half_width = 0
msms_min_feature_length = 1
msms_max_internal_gap = 0
msms_min_window_intensity = {t2}
msms_min_feature_intensity = 0
msms_iterations = 1
all_frames = false
"""


def peaks(d: str, ms_level: str) -> int:
    op = "=" if ms_level == "ms1" else "!="
    c = sqlite3.connect(d + "/analysis.tdf")
    v = c.execute(f"SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType{op}0").fetchone()[0]
    c.close()
    return int(v)


def read_t1(ms1_cfg: Path) -> int:
    m = re.search(r"(?m)^min_window_intensity\s*=\s*(\d+)", ms1_cfg.read_text())
    if not m:
        raise SystemExit(f"could not find min_window_intensity in {ms1_cfg}")
    return int(m.group(1))


def _kept_ms2(rep1: str, T: int, r2: int, tmp: str) -> float:
    out = f"{tmp}/c.d"
    shutil.rmtree(out, ignore_errors=True)
    subprocess.run(
        [str(DNOISE), rep1, out, "--ms1-polygon", "--denoise-msms",
         "--msms-mz-half-width", "0", "--msms-min-feature-length", "1",
         "--msms-max-internal-gap", "0", "--msms-iterations", "1",
         "--msms-min-window-intensity", str(T)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    k = peaks(out, "ms2") / r2
    shutil.rmtree(out, ignore_errors=True)
    return k


def calibrate_t2(rep1: str, target: float, r2: int, tmp: str, tol: float = 0.01) -> int:
    pts = []
    for T in LADDER:
        k = _kept_ms2(rep1, T, r2, tmp)
        pts.append((T, k))
        print(f"    T2={T:>5d}: keeps {k * 100:5.1f}% of MS2")
    pts.sort()
    lo = hi = None
    for (t0, k0), (t1, k1) in zip(pts, pts[1:]):
        if k0 >= target >= k1:
            lo, hi = t0, t1
            break
    best = min(pts, key=lambda p: abs(p[1] - target))
    if lo is None:
        return best[0]
    while hi - lo > 1:
        mid = (lo + hi) // 2
        k = _kept_ms2(rep1, mid, r2, tmp)
        print(f"    bisect T2={mid:>5d}: keeps {k * 100:5.1f}% (target {target * 100:.1f}%)")
        if abs(k - target) < abs(best[1] - target):
            best = (mid, k)
        if k > target:
            lo = mid
        else:
            hi = mid
        if abs(k - target) <= tol:
            break
    return best[0]


def main() -> int:
    if not DNOISE.exists():
        raise SystemExit(f"dnoise binary missing: {DNOISE} (cargo build --release)")
    for ds, cfg in GRADS.items():
        rep1 = sorted(glob.glob(f"{ROOT}/data/{ds}/raw/*REP1.d"))[0]
        r2 = peaks(rep1, "ms2")
        t1 = read_t1(cfg["ms1_cfg"])
        print(f"\n=== {ds}: T1={t1} (reused from {cfg['ms1_cfg'].name}) ===")

        with tempfile.TemporaryDirectory() as tmp:
            # streak MS1+MS/MS reference on REP1 (matches 04_denoise.sh's msms arm)
            streak_out = f"{tmp}/streak.d"
            subprocess.run(
                [str(DNOISE), rep1, streak_out, "--config", str(STREAK_CFG),
                 "--denoise-msms", "--ms1-polygon"],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
            target = peaks(streak_out, "ms2") / r2
            shutil.rmtree(streak_out, ignore_errors=True)
            print(f"  streak msms arm keeps {target * 100:.1f}% of MS2 (target)")

            print("  calibrating T2 (MS2):")
            t2 = calibrate_t2(rep1, target, r2, tmp)
        print(f"  -> T1={t1} (reused), T2={t2}")

        suffix = ".15min" if ds == "dda_15min" else ""
        cfg["out_cfg"].write_text(TEMPLATE.format(suffix=suffix, t1=t1, t2=t2))
        print(f"  wrote {cfg['out_cfg'].name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
