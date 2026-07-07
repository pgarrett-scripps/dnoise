#!/usr/bin/env python3
"""Calibrate matched per-point intensity thresholds T1 (MS1) and T2 (MS1+MS2)
for the main HYE diaPASEF benchmark (dia_5min / dia_15min), the diaPASEF analog
of 55_recalib_intensity_msms.py and the main-dataset port of
40_dia_ups_intensity_calib.py.

Unlike the UPS datasets, data/dia_5min and data/dia_15min normally hold only
raw/ on disk (the streak denoised/denoised_msms arms are built, searched, and
then deleted to save space -- see 04_denoise.sh / 11_diann.sh). So this script
rebuilds temporary single-REP1 streak references on the fly (--dia-ms1-window
for the MS1 target, --dia-ms1-window --dia-window --denoise-msms for the MS2
target), measures their kept-fraction, and deletes them immediately -- it never
holds more than one single-file .d copy on disk at a time.

For each DIA gradient:
  T1 (MS1)     -> ladder + bisection (--dia-ms1-window --min-window-intensity T)
                  matched to the streak MS1 arm's kept MS1 fraction on REP1.
  T2 (MS1+MS2) -> ladder + bisection (--dia-ms1-window --dia-window
                  --denoise-msms --msms-min-window-intensity T, degenerate MS2
                  params) matched to the streak msms arm's kept MS2 fraction.
T1 and T2 are calibrated independently (MS1/MS2 point counts don't interact;
same reasoning as 36_ups_intensity_calib.py / 40_dia_ups_intensity_calib.py).

Writes:
  config/dnoise.intensity.dia_5min.toml        (MS1-only, T1)
  config/dnoise.intensity.dia_15min.toml
  config/dnoise.intensity_msms.dia_5min.toml   (MS1+MS2, T1+T2)
  config/dnoise.intensity_msms.dia_15min.toml

Usage: uv run scripts/56_recalib_intensity_dia.py
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
STREAK_CFG = ROOT / "config" / "dnoise.toml"
DATASETS = ["dia_5min", "dia_15min"]
LADDER = [20, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800]

MS1_TEMPLATE = """\
# Strict intensity-threshold baseline, MS1-only diaPASEF (control arm for the
# streak `denoised` arm). Collapses the vertical filter to a pure per-point
# threshold via mz_half_width=0, min_feature_length=1, max_internal_gap=0,
# iterations=1, halo=false (see dnoise.intensity.toml for the DDA rationale;
# identical trick, DIA gate is --dia-ms1-window instead of --ms1-polygon).
#
# T = {t1} was calibrated by 56_recalib_intensity_dia.py to match the streak
# MS1 arm's (--dia-ms1-window) MS1-point removal on this gradient's REP1.

mz_half_width = 0
min_feature_length = 1
max_internal_gap = 0
min_window_intensity = {t1}
min_feature_intensity = 0
iterations = 1
halo = false
all_frames = false
"""

MSMS_TEMPLATE = """\
# Strict intensity-threshold baseline, MS1+MS2 diaPASEF (control arm for the
# streak `msms` arm). Both frame types collapsed to a pure per-point threshold
# (see dnoise.intensity.dia_5min.toml / dnoise.intensity_msms.toml for the
# rationale). DIA gates: --dia-ms1-window (MS1) + --dia-window (MS2).
#
# T1 = {t1}, T2 = {t2} were calibrated by 56_recalib_intensity_dia.py to match
# the streak `denoised`/`msms` arms' MS1/MS2-point removal on this gradient's
# REP1 (independently: MS1 and MS2 point counts don't interact).

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


def _kept(rep1: str, level: str, flags: list[str], raw_peaks: int, tmp: str) -> float:
    out = f"{tmp}/c.d"
    shutil.rmtree(out, ignore_errors=True)
    subprocess.run([str(DNOISE), rep1, out, *flags],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    k = peaks(out, level) / raw_peaks
    shutil.rmtree(out, ignore_errors=True)
    return k


def calibrate(rep1: str, level: str, flags_for_T, target: float, raw_peaks: int,
              tmp: str, tol: float = 0.01) -> int:
    pts = []
    for T in LADDER:
        k = _kept(rep1, level, flags_for_T(T), raw_peaks, tmp)
        pts.append((T, k))
        print(f"    T={T:>6d}: keeps {k * 100:5.1f}% of {level.upper()}")
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
        k = _kept(rep1, level, flags_for_T(mid), raw_peaks, tmp)
        print(f"    bisect T={mid:>6d}: keeps {k * 100:5.1f}% (target {target * 100:.1f}%)")
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
    for ds in DATASETS:
        rep1 = sorted(glob.glob(f"{ROOT}/data/{ds}/raw/*REP1.d"))[0]
        r1, r2 = peaks(rep1, "ms1"), peaks(rep1, "ms2")
        print(f"\n===== {ds} (REP1: {Path(rep1).name}) =====")

        with tempfile.TemporaryDirectory() as tmp:
            # --- streak MS1 target (--dia-ms1-window, matches 04_denoise.sh ms1 arm) ---
            streak_ms1 = f"{tmp}/streak_ms1.d"
            subprocess.run([str(DNOISE), rep1, streak_ms1, "--config", str(STREAK_CFG),
                             "--dia-ms1-window"],
                            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            tgt1 = peaks(streak_ms1, "ms1") / r1
            shutil.rmtree(streak_ms1, ignore_errors=True)
            print(f"  streak ms1 arm keeps {tgt1 * 100:.1f}% of MS1 (target T1)")

            # --- streak MS1+MS2 target (--dia-ms1-window --dia-window --denoise-msms,
            #     matches 04_denoise.sh msms arm) ---
            streak_msms = f"{tmp}/streak_msms.d"
            subprocess.run([str(DNOISE), rep1, streak_msms, "--config", str(STREAK_CFG),
                             "--dia-ms1-window", "--dia-window", "--denoise-msms"],
                            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            tgt2 = peaks(streak_msms, "ms2") / r2
            shutil.rmtree(streak_msms, ignore_errors=True)
            print(f"  streak msms arm keeps {tgt2 * 100:.1f}% of MS2 (target T2)")

            print("  calibrating T1 (MS1, --dia-ms1-window):")
            t1 = calibrate(
                rep1, "ms1",
                lambda T: ["--dia-ms1-window", "--min-window-intensity", str(T)],
                tgt1, r1, tmp)
            print("  calibrating T2 (MS2, --dia-ms1-window --dia-window):")
            t2 = calibrate(
                rep1, "ms2",
                lambda T: ["--dia-ms1-window", "--dia-window", "--denoise-msms",
                           "--msms-mz-half-width", "0", "--msms-min-feature-length", "1",
                           "--msms-max-internal-gap", "0", "--msms-iterations", "1",
                           "--msms-min-window-intensity", str(T)],
                tgt2, r2, tmp)
        print(f"  -> T1={t1} (MS1), T2={t2} (MS2)")

        (ROOT / "config" / f"dnoise.intensity.{ds}.toml").write_text(MS1_TEMPLATE.format(t1=t1))
        (ROOT / "config" / f"dnoise.intensity_msms.{ds}.toml").write_text(MSMS_TEMPLATE.format(t1=t1, t2=t2))
        print(f"  wrote dnoise.intensity.{ds}.toml + dnoise.intensity_msms.{ds}.toml")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
