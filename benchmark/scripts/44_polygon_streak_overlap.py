#!/usr/bin/env python3
"""Decompose ddaPASEF MS1 reduction into the selection-polygon gate vs. the
streak (vertical+halo) filter, to show they are NOT additive contributions.

Motivation (reviewer rebuttal): the headline "~83-85% of MS1 removed" bundles the
polygon gate and the structural filter. A reviewer reads this as "half the number
is a trivial geometric crop." But the gate and the streak filter overlap: the
streak filter run on the FULL frame (no gate) already removes most of the
out-of-polygon region, because that region (the singly-charged trend line the
PASEF method excludes from precursor selection) is mostly scattered, non-streaky
signal. The gate is therefore largely a fast shortcut for what the streak prior
removes anyway; its only UNIQUE contribution is the real, vertical-streak ion
signal that happens to fall outside the selection polygon (never-fragmented
single-charge ions).

We measure four MS1 point counts per run, all from the SAME binary/config so they
are comparable, then derive the Venn decomposition by inclusion-exclusion:

  T       total MS1 points              (raw)
  R_both  removed by polygon + streak   (the default DDA arm: --ms1-polygon)
  R_strk  removed by streak alone       (no --ms1-polygon)
  R_poly  removed by polygon alone      (--ms1-polygon, streak neutralized:
                                         --min-feature-length 1 --no-halo)

Since applying both removes the UNION of the two removal sets:
  overlap     = R_poly + R_strk - R_both     (out-of-polygon ALSO killed by streak)
  gate_unique = R_both - R_strk  (= R_poly - overlap)   out-of-polygon streak KEEPS
  strk_inpoly = R_both - R_poly  (= R_strk - overlap)   in-polygon noise streak kills
  kept        = T - R_both
These four partition T exactly.

Disk-light & idempotent: for each raw .d we denoise the three arms into a temp
dir, read SUM(NumPeaks) for MsMsType=0, then delete the temp .d. Per-file counts
are cached in results/_overlap/polygon_streak_overlap.csv so reruns skip done
files. Run:  uv run scripts/44_polygon_streak_overlap.py
"""

from __future__ import annotations

import shutil
import sqlite3
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DNOISE = ROOT.parent / "target" / "release" / "dnoise"
CONFIG = ROOT / "config" / "dnoise.toml"
DATASETS = ["dda_5min", "dda_15min"]
OUT = ROOT / "results" / "_overlap"
CSV = OUT / "polygon_streak_overlap.csv"
FIG = ROOT.parent / "paper" / "figures" / "si_polygon_streak_overlap.png"

# Colorblind-safe (Wong, Nat. Methods 2011).
C_KEPT = "#999999"   # signal kept by the full pipeline
C_INPOLY = "#0072B2"  # in-polygon noise removed by the streak filter
C_OVERLAP = "#56B4E9"  # out-of-polygon removed by streak too (gate/streak overlap)
C_GATE = "#E69F00"   # out-of-polygon kept by streak, removed only by the gate


def ms1_peaks(d: Path) -> int:
    c = sqlite3.connect(d / "analysis.tdf")
    n = c.execute(
        "SELECT COALESCE(SUM(NumPeaks),0) FROM Frames WHERE MsMsType=0"
    ).fetchone()[0]
    c.close()
    return int(n)


def denoise_count(raw: Path, extra: list[str]) -> int:
    """Denoise one raw .d with the given extra flags into a temp dir, return the
    MS1 peak count, and delete the temp output."""
    with tempfile.TemporaryDirectory(prefix="dnoise_overlap_") as tmp:
        out = Path(tmp) / raw.name
        subprocess.run(
            [str(DNOISE), str(raw), str(out), "--config", str(CONFIG), *extra],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return ms1_peaks(out)


def short(fname: str) -> str:
    import re

    m = re.search(r"Condition_([ABC])_REP(\d)", fname)
    return f"{m.group(1)}{m.group(2)}" if m else fname


def collect() -> pd.DataFrame:
    OUT.mkdir(parents=True, exist_ok=True)
    done = set()
    if CSV.exists():
        prev = pd.read_csv(CSV)
        done = set(zip(prev["dataset"], prev["file"]))
    else:
        prev = pd.DataFrame()

    rows = []
    for ds in DATASETS:
        raw_dir = ROOT / "data" / ds / "raw"
        for raw in sorted(raw_dir.glob("*.d")):
            key = (ds, raw.name)
            if key in done:
                continue
            print(f"[{ds}] {short(raw.name)} ...", flush=True)
            total = ms1_peaks(raw)
            both = total - denoise_count(raw, ["--ms1-polygon"])
            strk = total - denoise_count(raw, [])
            poly = total - denoise_count(
                raw, ["--ms1-polygon", "--min-feature-length", "1", "--no-halo"]
            )
            rows.append({
                "dataset": ds, "file": raw.name, "short": short(raw.name),
                "total": total, "R_both": both, "R_strk": strk, "R_poly": poly,
            })
            # Persist incrementally so the run is restartable.
            pd.concat([prev, pd.DataFrame(rows)], ignore_index=True).to_csv(
                CSV, index=False
            )
    return pd.read_csv(CSV)


def decompose(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["overlap"] = d["R_poly"] + d["R_strk"] - d["R_both"]
    d["gate_unique"] = d["R_both"] - d["R_strk"]
    d["strk_inpoly"] = d["R_both"] - d["R_poly"]
    d["kept"] = d["total"] - d["R_both"]
    return d


def main() -> int:
    if not DNOISE.exists():
        raise SystemExit(f"dnoise binary not found: {DNOISE} (cargo build --release)")
    df = decompose(collect())

    # Aggregate per gradient: sum point counts across runs, express as % of total.
    fig, axes = plt.subplots(1, len(DATASETS), figsize=(10, 4.2), sharey=True)
    if len(DATASETS) == 1:
        axes = [axes]
    summary = {}
    for ax, ds in zip(axes, DATASETS):
        g = df[df["dataset"] == ds]
        T = g["total"].sum()
        kept = g["kept"].sum() / T * 100
        inpoly = g["strk_inpoly"].sum() / T * 100
        overlap = g["overlap"].sum() / T * 100
        gate = g["gate_unique"].sum() / T * 100
        # streak-alone and gate-alone removal as % of all MS1
        strk_all = g["R_strk"].sum() / T * 100
        poly_all = g["R_poly"].sum() / T * 100
        both_all = g["R_both"].sum() / T * 100
        # fraction of the gate's removed points that the streak filter ALSO removes
        overlap_of_gate = g["overlap"].sum() / g["R_poly"].sum() * 100
        summary[ds] = {
            "n_runs": int(len(g)), "total_ms1": int(T),
            "pct_kept": kept, "pct_streak_inpoly": inpoly,
            "pct_overlap": overlap, "pct_gate_unique": gate,
            "pct_removed_streak_alone": strk_all,
            "pct_removed_polygon_alone": poly_all,
            "pct_removed_both": both_all,
            "pct_of_gate_region_also_removed_by_streak": overlap_of_gate,
        }

        # One stacked bar: kept | in-polygon streak | overlap | gate-only.
        segs = [
            ("kept (signal)", kept, C_KEPT),
            ("in-polygon noise\nremoved by streak", inpoly, C_INPOLY),
            ("out-of-polygon removed\nby streak too (overlap)", overlap, C_OVERLAP),
            ("out-of-polygon removed\nonly by gate", gate, C_GATE),
        ]
        left = 0.0
        for label, val, color in segs:
            ax.barh(0, val, left=left, color=color, edgecolor="white",
                    label=label if ax is axes[0] else None)
            if val > 3:
                ax.text(left + val / 2, 0, f"{val:.0f}%", ha="center", va="center",
                        fontsize=10, color="white", fontweight="bold")
            left += val
        grad = ds.replace("dda_", "").replace("min", "-min")
        ax.set_title(
            f"ddaPASEF {grad}\nstreak filter alone: {strk_all:.0f}%   "
            f"gate-unique: +{gate:.0f}%",
            fontsize=11,
        )
        ax.set_xlim(0, 100)
        ax.set_yticks([])
        ax.set_xlabel("% of total MS1 points")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(
        "MS1 reduction decomposed: the polygon gate and the streak filter overlap",
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=200, bbox_inches="tight")
    print(f"\nwrote {FIG}")

    # Print a compact summary for the SI prose.
    sm = pd.DataFrame(summary).T
    pd.set_option("display.float_format", lambda x: f"{x:.1f}")
    print("\n=== decomposition (% of total MS1, pooled across runs) ===")
    print(sm.to_string())
    (OUT / "polygon_streak_overlap_summary.csv").write_text(sm.to_csv())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
