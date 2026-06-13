#!/usr/bin/env python3
"""Quantify dnoise's data-volume reduction across three arms: raw, MS1-denoised
(data/denoised), and MS1+MS/MS-denoised (data/denoised_msms).

Writes a two-panel figure to results/analysis/:
  - per-file frame-binary size for the three arms (shows both stages + consistency)
  - aggregate peaks split into MS1 / MS/MS for the three arms (shows which stage
    removes what: MS1 denoising cuts the MS1 block; MS/MS denoising cuts the MS2
    block on top).
plus a CSV.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _metrics import stats

ROOT = Path(__file__).resolve().parents[1]
DATASET = os.environ.get("DATASET", "dda_5min")
DATA = ROOT / "data" / DATASET
RAW = DATA / "raw"
ARMS = {
    "raw": RAW,
    "MS1": DATA / "denoised",
    "MS1+MS/MS": DATA / "denoised_msms",
}
OUT = ROOT / "results" / DATASET / "analysis"
COLOR = {"raw": "#1f77b4", "MS1": "#d62728", "MS1+MS/MS": "#2ca02c"}


def short_name(fname: str) -> str:
    m = re.search(r"Condition_([ABC])_REP(\d)", fname)
    return f"{m.group(1)}{m.group(2)}" if m else fname


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for raw_d in sorted(RAW.glob("*.d")):
        rec = {"file": short_name(raw_d.name)}
        ok = True
        for arm, base in ARMS.items():
            d = base / raw_d.name
            if not d.is_dir():
                ok = False
                break
            ms1, ms2, b = stats(d)
            rec[f"ms1_{arm}"], rec[f"ms2_{arm}"], rec[f"bytes_{arm}"] = ms1, ms2, b
        if ok:
            rows.append(rec)
    df = pd.DataFrame(rows).sort_values("file").reset_index(drop=True)
    if df.empty:
        print("no matching .d across all three arms")
        return 1
    df.to_csv(OUT / "data_reduction.csv", index=False)

    arms = list(ARMS)
    # Aggregate totals per arm (GB, billions of peaks).
    agg = {a: {
        "bytes": df[f"bytes_{a}"].sum(),
        "ms1": df[f"ms1_{a}"].sum(),
        "ms2": df[f"ms2_{a}"].sum(),
    } for a in arms}
    raw_b = agg["raw"]["bytes"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Panel A: per-file frame-binary size, three bars per file.
    ax = axes[0]
    x = np.arange(len(df))
    w = 0.27
    for i, a in enumerate(arms):
        ax.bar(x + (i - 1) * w, df[f"bytes_{a}"] / 1e6, w, label=a, color=COLOR[a])
    ax.set_xticks(x)
    ax.set_xticklabels(df["file"], rotation=90, fontsize=7)
    ax.set_ylabel("analysis.tdf_bin (MB)")
    pct = {a: 100 * (1 - agg[a]["bytes"] / raw_b) for a in arms}
    ax.set_title(f"Frame binary size per run\n(overall −{pct['MS1']:.0f}% MS1, −{pct['MS1+MS/MS']:.0f}% MS1+MS/MS)")
    ax.legend()

    # Panel B: aggregate peaks, stacked MS1 (solid) + MS/MS (hatched), per arm.
    ax = axes[1]
    xb = np.arange(len(arms))
    ms1 = [agg[a]["ms1"] / 1e9 for a in arms]
    ms2 = [agg[a]["ms2"] / 1e9 for a in arms]
    ax.bar(xb, ms1, 0.6, label="MS1 peaks", color="#4c72b0")
    ax.bar(xb, ms2, 0.6, bottom=ms1, label="MS/MS peaks", color="#dd8452")
    for j, a in enumerate(arms):
        tot = ms1[j] + ms2[j]
        ax.text(j, tot + 0.08, f"{tot:.2f}B\n(−{100*(1-tot/((agg['raw']['ms1']+agg['raw']['ms2'])/1e9)):.0f}%)",
                ha="center", va="bottom", fontsize=8)
    ax.set_xticks(xb)
    ax.set_xticklabels(arms)
    ax.set_ylabel("peaks (billions)")
    ax.set_ylim(0, (agg["raw"]["ms1"] + agg["raw"]["ms2"]) / 1e9 * 1.15)
    ax.set_title("Peaks by level (MS1 denoising cuts MS1; MS/MS denoising cuts MS/MS)")
    ax.legend()

    fig.suptitle(f"dnoise data-volume reduction ({len(df)} runs, default settings)")
    fig.tight_layout()
    fig.savefig(OUT / "data_reduction.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("aggregate (raw -> MS1 -> MS1+MS/MS):")
    print(f"  binary GB : {agg['raw']['bytes']/1e9:.2f} -> {agg['MS1']['bytes']/1e9:.2f} -> {agg['MS1+MS/MS']['bytes']/1e9:.2f}")
    print(f"  MS1 peaks : {agg['raw']['ms1']/1e6:.0f}M -> {agg['MS1']['ms1']/1e6:.0f}M -> {agg['MS1+MS/MS']['ms1']/1e6:.0f}M")
    print(f"  MS2 peaks : {agg['raw']['ms2']/1e6:.0f}M -> {agg['MS1']['ms2']/1e6:.0f}M -> {agg['MS1+MS/MS']['ms2']/1e6:.0f}M")
    print(f"\nWrote {OUT}/data_reduction.png and data_reduction.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
