#!/usr/bin/env python3
"""Recalibrate the intensity-threshold control's T to match the streak arm's MS1
removal, WITH the polygon gate applied first (the streak/denoised arm now gates
on the ddaPASEF selection polygon, so the matched brightness-cut control must too).

For each DDA gradient: take REP1, apply --ms1-polygon + a ladder of per-point
thresholds T, measure total MS1 removal, and pick the T whose removal matches the
streak (denoised) arm's removal for that gradient. Writes T into the per-gradient
intensity config so 04_denoise.sh (intensity arm, which already adds --ms1-polygon)
reproduces a matched-reduction brightness cut.
"""

from __future__ import annotations

import glob
import re
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DNOISE = (ROOT / ".." / "target" / "release" / "dnoise").resolve()

# gradient -> (intensity config, REP1 glob, raw glob)
GRADS = {
    "dda_5min": ROOT / "config" / "dnoise.intensity.toml",
    "dda_15min": ROOT / "config" / "dnoise.intensity.15min.toml",
}
LADDER = [40, 80, 120, 180, 260, 380]


def ms1(d: str) -> int:
    c = sqlite3.connect(d + "/analysis.tdf")
    v = c.execute("SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType=0").fetchone()[0]
    c.close()
    return int(v)


def streak_keep_frac(ds: str) -> float:
    raw = sum(ms1(d) for d in glob.glob(f"{ROOT}/data/{ds}/raw/*.d"))
    den = sum(ms1(d) for d in glob.glob(f"{ROOT}/data/{ds}/denoised/*.d"))
    return den / raw


def main() -> int:
    for ds, cfg in GRADS.items():
        rep1 = sorted(glob.glob(f"{ROOT}/data/{ds}/raw/*REP1.d"))[0]
        raw1 = ms1(rep1)
        target = streak_keep_frac(ds)  # streak's kept fraction of raw MS1
        print(f"\n=== {ds}: streak keeps {target*100:.1f}% of MS1 (target) ===")
        pts = []
        for T in LADDER:
            out = f"/tmp/recal_{ds}.d"
            shutil.rmtree(out, ignore_errors=True)
            subprocess.run([str(DNOISE), rep1, out, "--config", str(cfg),
                            "--ms1-polygon", "--min-window-intensity", str(T)],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            keep = ms1(out) / raw1
            shutil.rmtree(out, ignore_errors=True)
            pts.append((T, keep))
            print(f"  T={T:>4d}: keeps {keep*100:.1f}% of MS1")

        # Pick the T whose kept-fraction is closest to the streak target; refine by
        # linear interpolation between the two bracketing ladder points.
        pts.sort()
        best_T = min(pts, key=lambda p: abs(p[1] - target))[0]
        for (t0, k0), (t1, k1) in zip(pts, pts[1:]):
            if (k0 - target) * (k1 - target) <= 0 and k0 != k1:
                best_T = round(t0 + (target - k0) * (t1 - t0) / (k1 - k0))
                break
        print(f"  -> chosen T = {best_T}")

        # Write T into the config's min_window_intensity line.
        text = cfg.read_text()
        new = re.sub(r"(?m)^min_window_intensity\s*=.*$",
                     f"min_window_intensity = {best_T}", text)
        cfg.write_text(new)
        print(f"  updated {cfg.name}: min_window_intensity = {best_T}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
