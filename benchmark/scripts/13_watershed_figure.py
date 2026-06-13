#!/usr/bin/env python3
"""Figure 6 — the watershed centroider arm.

A closing highlight figure: run the MS1 vertical+halo denoiser plus the watershed
centroider (config/dnoise.watershed.toml, MS/MS untouched), then contrast it
against the original raw data on three axes:

  - Compression: aggregate MS1 peaks and analysis.tdf_bin size. The centroider
    collapses each MS1 frame to a handful of intensity-weighted points, so this
    is the dramatic panel (MS1 reduced to a tiny fraction).
  - Identifications: PSM / peptide / protein at 1% FDR. MS/MS is untouched, so
    DDA IDs are expected to be essentially unchanged.
  - LFQ: observed median log2(A/B) per species vs the known HYE ratios, plus the
    median protein CV. This is where the aggressive MS1 reduction shows its cost.

Reads:
  data/<dataset>/{raw,watershed}/*.d         (compression)
  results/<dataset>/{original,watershed}/    (IDs, LFQ)
Writes:
  results/<dataset>/analysis/fig6_watershed.{png,csv}
  paper/figures/fig6_watershed.png           (if paper/figures exists)

The watershed arm is optional; if its data/results are absent this prints how to
generate them and exits 0 (so it never breaks the main pipeline).
"""

from __future__ import annotations

import os
import sqlite3
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

from _metrics import (
    EXPECTED,
    SPECIES,
    id_metrics,
    lfq_metrics,
    lfq_table,
    pair_accuracy,
    stats,
)

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
DATASET = os.environ.get("DATASET", "dda_5min")
DATA = ROOT / "data" / DATASET
RESULTS = ROOT / "results" / DATASET
OUT = RESULTS / "analysis"
PAPER_FIGS = ROOT.parent / "paper" / "figures"

RAW_D = DATA / "raw"
WSHED_D = DATA / "watershed"
RES_ORIG = RESULTS / "original"
RES_WSHED = RESULTS / "watershed"

# Frame point-cloud panels reuse the same frame/zoom as Figure 1 (08_frame_figures).
DUMP = REPO / "target" / "release" / "examples" / "dump_frame"
FRAME_FILE = "LFQ_Ultra2_PASEF_5min_50ng_Condition_A_REP1.d"
MZ_HALF = 5.0      # zoom is 10 m/z wide (matches Figure 1)
K0_HALF = 0.075    # zoom is 0.15 1/K0 wide (matches Figure 1)
CMAP = "viridis"

ARM_COLOR = {"original": "#1f77b4", "watershed": "#9467bd"}


def most_intense_ms1_index(d: Path) -> int:
    """0-based timsrust index of the MS1 frame with the largest summed intensity."""
    c = sqlite3.connect(d / "analysis.tdf")
    ids = [r[0] for r in c.execute("SELECT Id FROM Frames ORDER BY Id")]
    best = c.execute(
        "SELECT Id FROM Frames WHERE MsMsType=0 ORDER BY SummedIntensities DESC LIMIT 1"
    ).fetchone()[0]
    c.close()
    return ids.index(best)


def dump_frame(d: Path, idx: int, out: Path) -> pd.DataFrame:
    subprocess.run([str(DUMP), str(d), str(idx), str(out)], check=True, stderr=subprocess.DEVNULL)
    return pd.read_csv(out)


def frame_view() -> dict | None:
    """Dump the watershed-centroided version of Figure 1's frame, plus the raw
    frame's zoom centre / axis extents (so the panels are directly comparable to
    Figure 1). Returns None if dump_frame or the input .d folders are missing."""
    raw_d = RAW_D / FRAME_FILE
    wshed_d = WSHED_D / FRAME_FILE
    if not DUMP.is_file() or not raw_d.is_dir() or not wshed_d.is_dir():
        return None
    idx = most_intense_ms1_index(raw_d)
    with tempfile.TemporaryDirectory() as tmp:
        raw = dump_frame(raw_d, idx, Path(tmp) / "raw.csv")
        wshed = dump_frame(wshed_d, idx, Path(tmp) / "wshed.csv")
    peak = raw.loc[raw["intensity"].idxmax()]  # zoom centre = raw frame's brightest point
    return {
        "wshed": wshed,
        "n_raw": len(raw),
        "mz0": float(peak["mz"]),
        "k0_0": float(peak["one_over_k0"]),
        "vmax": float(raw["intensity"].max()),
        "mz_edges": np.linspace(raw["mz"].min(), raw["mz"].max(), 400),
        "k0_edges": np.linspace(raw["one_over_k0"].min(), raw["one_over_k0"].max(), 300),
    }


def frame_hist(ax, df, mz_edges, k0_edges, norm):
    h = ax.hist2d(df["mz"], df["one_over_k0"], bins=[mz_edges, k0_edges],
                  weights=df["intensity"], norm=norm, cmap=CMAP)
    ax.set_xlabel("m/z")
    ax.set_ylabel("ion mobility (1/K0)")
    return h[3]


def zoom_scatter(ax, df, mz0, k0_0, norm):
    sub = df[df["mz"].between(mz0 - MZ_HALF, mz0 + MZ_HALF)
             & df["one_over_k0"].between(k0_0 - K0_HALF, k0_0 + K0_HALF)]
    sc = ax.scatter(sub["mz"], sub["one_over_k0"], c=sub["intensity"],
                    norm=norm, cmap=CMAP, s=12, edgecolors="none")
    ax.set_xlim(mz0 - MZ_HALF, mz0 + MZ_HALF)
    ax.set_ylim(k0_0 - K0_HALF, k0_0 + K0_HALF)
    ax.set_xlabel("m/z")
    ax.set_ylabel("ion mobility (1/K0)")
    return sc, sub


def aggregate_compression() -> dict | None:
    """Sum MS1 peaks / MS2 peaks / tdf_bin bytes over .d folders present in both
    the raw and watershed arms."""
    if not RAW_D.is_dir() or not WSHED_D.is_dir():
        return None
    agg = {a: {"ms1": 0, "ms2": 0, "bytes": 0} for a in ("original", "watershed")}
    n = 0
    for raw_d in sorted(RAW_D.glob("*.d")):
        w_d = WSHED_D / raw_d.name
        if not w_d.is_dir():
            continue
        for arm, d in (("original", raw_d), ("watershed", w_d)):
            ms1, ms2, b = stats(d)
            agg[arm]["ms1"] += ms1
            agg[arm]["ms2"] += ms2
            agg[arm]["bytes"] += b
        n += 1
    if n == 0:
        return None
    agg["_n"] = n
    return agg


def plot(agg: dict, lfq: dict, prot_orig, prot_wshed, acc: list[dict], fv: dict | None) -> None:
    from matplotlib.patches import Patch  # local import keeps the panels self-contained
    w = 0.38

    # Two rows when frame panels are available (top: watershed frame point clouds,
    # bottom: the four quantitative panels); otherwise just the bottom row.
    if fv is not None:
        fig = plt.figure(figsize=(18, 9.8))
        gs = fig.add_gridspec(2, 4, height_ratios=[1.15, 1.0], hspace=0.32, wspace=0.32)
        ax_full = fig.add_subplot(gs[0, 0:2])
        ax_zoom = fig.add_subplot(gs[0, 2:4])
        axes = [fig.add_subplot(gs[1, j]) for j in range(4)]

        # ---- Top panels: the SAME frame as Figure 1, watershed-centroided only ----
        wdf = fv["wshed"]
        full_norm = LogNorm(vmin=1, vmax=fv["vmax"])
        im = frame_hist(ax_full, wdf, fv["mz_edges"], fv["k0_edges"], full_norm)
        ax_full.add_patch(Rectangle(
            (fv["mz0"] - MZ_HALF, fv["k0_0"] - K0_HALF), 2 * MZ_HALF, 2 * K0_HALF,
            fill=False, edgecolor="red", lw=1.2))
        ax_full.set_title(f"watershed centroids — full MS1 frame ({len(wdf):,} points; "
                          f"raw had {fv['n_raw']:,})")
        fig.colorbar(im, ax=ax_full, label="summed intensity", shrink=0.85)
        sc, sub = zoom_scatter(ax_zoom, wdf, fv["mz0"], fv["k0_0"], full_norm)
        ax_zoom.set_title(f"watershed centroids — zoom ({len(sub):,} points)")
        fig.colorbar(sc, ax=ax_zoom, label="intensity", shrink=0.85)
    else:
        fig, axes = plt.subplots(1, 4, figsize=(19, 4.8))

    # ---- Panel A: compression (MS1 peaks + binary size, % of raw retained) ----
    ax = axes[0]
    raw_ms1 = agg["original"]["ms1"]
    w_ms1 = agg["watershed"]["ms1"]
    raw_b = agg["original"]["bytes"]
    w_b = agg["watershed"]["bytes"]
    metrics = ["MS1 peaks", "tdf_bin size"]
    orig_pct = [100.0, 100.0]
    wshed_pct = [100 * w_ms1 / raw_ms1, 100 * w_b / raw_b]
    x = np.arange(len(metrics))
    ax.bar(x - w / 2, orig_pct, w, label="original", color=ARM_COLOR["original"])
    bars = ax.bar(x + w / 2, wshed_pct, w, label="watershed", color=ARM_COLOR["watershed"])
    ax.bar_label(bars, fmt="%.1f%%", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylabel("% of raw retained")
    ax.set_ylim(0, 115)
    ax.set_title(
        f"Compression ({agg['_n']} runs)\n"
        f"MS1 {raw_ms1/1e9:.2f}B → {w_ms1/1e6:.0f}M peaks · "
        f"{raw_b/1e9:.1f} → {w_b/1e9:.1f} GB"
    )
    ax.legend()

    # ---- Panel B: proteins quantified (the only ID-level count expected to move) ----
    ax = axes[1]
    nq = [lfq["original"].get("n_quantified", 0), lfq["watershed"].get("n_quantified", 0)]
    bars = ax.bar([0, 1], nq, 0.6, color=[ARM_COLOR["original"], ARM_COLOR["watershed"]])
    ax.bar_label(bars, fmt="%d", fontsize=9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["original", "watershed"], rotation=20)
    ax.set_ylabel("proteins quantified @ 1% FDR")
    ax.set_ylim(0, max(nq) * 1.18 if max(nq) else 1)
    ax.set_title("Proteins quantified\n(IDs otherwise unchanged)")

    # ---- Panel C: LFQ accuracy across the dynamic range (observed vs expected) ----
    ax = axes[2]
    lim = (-3.5, 4.0)
    ax.plot(lim, lim, color="gray", ls="--", lw=1, zorder=0, label="ideal")
    marker = {"HUMAN": "o", "YEAST": "s", "ECOLI": "^"}
    acc_df = pd.DataFrame(acc)
    for arm in ("original", "watershed"):
        d = acc_df[acc_df["arm"] == arm] if not acc_df.empty else acc_df
        for sp in SPECIES:
            ds = d[d["species"] == sp] if not d.empty else d
            if ds.empty:
                continue
            ax.scatter(ds["expected"], ds["observed"], color=ARM_COLOR[arm],
                       marker=marker[sp], s=45, alpha=0.85,
                       label=arm if sp == "HUMAN" else None)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("expected log2 ratio")
    ax.set_ylabel("observed median log2 ratio")
    ax.set_title("LFQ accuracy across dynamic range\n(A/B, A/C, B/C; ○ human △ ecoli □ yeast)")
    ax.legend(fontsize=8)

    # ---- Panel D: LFQ ratio distributions (violins, A/B per species) ----
    ax = axes[3]
    arms = [("original", prot_orig), ("watershed", prot_wshed)]
    vw = 0.34  # violin half-slot width per arm
    offsets = {"original": -vw / 2, "watershed": vw / 2}
    for arm, prot in arms:
        if prot is None or prot.empty:
            continue
        data, pos = [], []
        for i, sp in enumerate(SPECIES):
            vals = prot[prot["species"] == sp]["log2_ratio"].dropna().values
            if len(vals) >= 2:  # violinplot needs a distribution
                data.append(vals)
                pos.append(i + offsets[arm])
        if not data:
            continue
        vp = ax.violinplot(data, positions=pos, widths=vw * 0.9,
                           showmedians=True, showextrema=False)
        for body in vp["bodies"]:
            body.set_facecolor(ARM_COLOR[arm])
            body.set_alpha(0.6)
            body.set_edgecolor("black")
            body.set_linewidth(0.5)
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(1.2)
    for i, sp in enumerate(SPECIES):
        ax.hlines(EXPECTED[sp], i - 0.45, i + 0.45, color="gray", lw=2, ls="--", zorder=5)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks(np.arange(len(SPECIES)))
    ax.set_xticklabels(SPECIES)
    ax.set_ylabel("log2(A/B) protein ratio")
    ax.set_ylim(-5, 4)
    cv_o = lfq["original"].get("median_cv", float("nan"))
    cv_w = lfq["watershed"].get("median_cv", float("nan"))
    ax.set_title(f"LFQ ratio distributions (dashed = expected)\nmedian CV: orig {cv_o:.3f} · wshed {cv_w:.3f}")
    ax.legend(handles=[Patch(facecolor=ARM_COLOR[a], alpha=0.6, edgecolor="black", label=a)
                       for a, _ in arms] + [plt.Line2D([], [], color="gray", ls="--", label="expected")])

    fig.suptitle("Watershed centroider: aggressive MS1 compression, IDs preserved, LFQ trade-off",
                 fontsize=13)
    if fv is None:
        fig.tight_layout()  # gridspec layout already sets its own spacing
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig6_watershed.png", dpi=150, bbox_inches="tight")
    if PAPER_FIGS.is_dir():
        fig.savefig(PAPER_FIGS / "fig6_watershed.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    missing = []
    if not WSHED_D.is_dir() or not any(WSHED_D.glob("*.d")):
        missing.append(f"watershed .d folders ({WSHED_D})")
    if not (RES_WSHED / "results.sage.tsv").is_file():
        missing.append(f"watershed Sage results ({RES_WSHED})")
    if missing:
        print("Figure 6 needs the watershed arm, which is not present yet:")
        for m in missing:
            print(f"  - missing: {m}")
        print("\nGenerate it with:")
        print(f"  DATASET={DATASET} just watershed   # denoise (wshed arm) + Sage")
        print(f"  DATASET={DATASET} just fig6        # then this figure")
        return 0

    agg = aggregate_compression()
    if agg is None:
        print("no .d folders shared between raw and watershed arms")
        return 1

    prot_orig = lfq_table(RES_ORIG)
    prot_wshed = lfq_table(RES_WSHED)
    lfq = {"original": {**id_metrics(RES_ORIG), **lfq_metrics(prot_orig)},
           "watershed": {**id_metrics(RES_WSHED), **lfq_metrics(prot_wshed)}}
    # Observed-vs-expected ratios across all pairs/species, for the accuracy panel.
    acc = pair_accuracy("original", RES_ORIG) + pair_accuracy("watershed", RES_WSHED)

    # Watershed-centroided point clouds of Figure 1's frame (top panels).
    fv = frame_view()
    if fv is None:
        print("note: dump_frame binary or input .d missing — skipping frame panels\n"
              "  build with: cargo build --release --example dump_frame")

    plot(agg, lfq, prot_orig, prot_wshed, acc, fv)

    # CSV summary alongside the figure.
    rows = []
    for arm in ("original", "watershed"):
        rows.append({
            "arm": arm,
            "ms1_peaks": agg[arm]["ms1"],
            "ms2_peaks": agg[arm]["ms2"],
            "tdf_bin_bytes": agg[arm]["bytes"],
            **{k: lfq[arm].get(k) for k in
               ("n_psm", "n_peptide", "n_protein", "n_quantified", "median_cv")},
            **{f"median_log2_{sp}": lfq[arm].get(f"median_log2_{sp}") for sp in SPECIES},
        })
    df = pd.DataFrame(rows).set_index("arm")
    df.to_csv(OUT / "fig6_watershed.csv")

    pd.set_option("display.width", 200, "display.max_columns", 100)
    print(df.T)
    print(f"\nMS1 peaks retained: {100*agg['watershed']['ms1']/agg['original']['ms1']:.2f}%  "
          f"({agg['original']['ms1']/1e9:.2f}B → {agg['watershed']['ms1']/1e6:.0f}M)")
    print(f"Wrote {OUT}/fig6_watershed.png and fig6_watershed.csv")
    if PAPER_FIGS.is_dir():
        print(f"Wrote {PAPER_FIGS}/fig6_watershed.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
