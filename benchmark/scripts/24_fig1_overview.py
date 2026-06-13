#!/usr/bin/env python3
"""Fig 1 — algorithm overview schematic (draft).

Left: the structural prior — in (m/z x ion-mobility) space a real ion is a
vertical streak (one m/z, spread only over its mobility peak) while noise is
scattered; the vertical filter keeps points in long, intense runs within a small
TOF-index window.
Right: the dnoise pipeline as a flow of stages (core = vertical + halo on MS1;
optional = MS/MS denoise and the two centroiders), writing a native type-2 .d.

This is a programmatic draft; the final figure can be redrawn in a vector tool.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

FIGS = Path(__file__).resolve().parents[1].parent / "paper" / "figures"


def concept(ax) -> None:
    rng = np.random.default_rng(0)
    # Physical-unit view: a narrow m/z slice vs. ion mobility. A real ion sits at
    # one m/z and spreads only over its mobility peak (a vertical streak); noise is
    # scattered. The filter keeps points in long, intense runs within a tight m/z
    # (TOF-index) window.
    mz_lo, mz_hi = 500.0, 501.0
    im_lo, im_hi = 0.75, 0.85
    mz0 = 500.5  # the real ion's m/z

    # scattered noise across the panel
    ax.scatter(rng.uniform(mz_lo, mz_hi, 120), rng.uniform(im_lo, im_hi, 120),
               s=10, c="#bbbbbb", label="noise (scattered)")
    # a couple of short noise runs (a few points sharing an m/z) — too short to keep
    for cx in (500.22, 500.82):
        ax.scatter([cx] * 3, rng.uniform(im_lo + 0.012, im_hi - 0.012, 3), s=12, c="#999999")
    # a real ion: vertical streak at mz0 spanning its mobility peak (uniform marker size)
    streak_y = np.linspace(0.782, 0.818, 12)
    ax.scatter([mz0] * len(streak_y), streak_y, s=40,
               c="#1f77b4", label="real ion (vertical streak)", zorder=3)
    # filter window around the streak's m/z column (TOF-index half-width; widened
    # here for visibility). Drawn just inside the panel so its label has headroom.
    half = 0.06
    ax.add_patch(Rectangle((mz0 - half, im_lo + 0.002), 2 * half, (im_hi - im_lo) - 0.004,
                           fill=False, edgecolor="#d62728", lw=1.6, ls="--"))
    ax.set_xlim(mz_lo, mz_hi); ax.set_ylim(im_lo, im_hi + 0.020)
    # window label sits in the top headroom, clear of the plotted points
    ax.text(mz0, im_hi + 0.004, "TOF window (±mz_half_width)",
            color="#d62728", fontsize=8, ha="center", va="bottom")
    # annotate the kept streak from empty space to the right
    ax.annotate("kept: long,\nintense run", xy=(mz0 + half, 0.806), xytext=(500.70, 0.828),
                fontsize=8, color="#1f77b4", ha="left",
                arrowprops=dict(arrowstyle="->", color="#1f77b4"))
    ax.set_xlabel("m/z"); ax.set_ylabel("ion mobility (1/K0)")
    ax.set_title("Structural prior", fontsize=11)
    ax.legend(loc="lower left", fontsize=7, framealpha=0.9)


def box(ax, x, y, w, h, text, fc, fontsize=9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.04",
                                fc=fc, ec="black", lw=1.0))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=14, lw=1.2, color="black"))


def pipeline(ax) -> None:
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.set_title("dnoise pipeline (MS1; MS/MS optional)", fontsize=11)
    # core path
    box(ax, 0.2, 6.2, 1.9, 1.4, "Bruker .d\n(type-2)", "#eaeaea")
    box(ax, 2.6, 6.2, 2.1, 1.4, "vertical-IM\nfeature filter", "#cfe2f3")
    box(ax, 5.1, 6.2, 2.0, 1.4, "horizontal-\nhalo filter", "#cfe2f3")
    box(ax, 7.5, 6.2, 2.2, 1.4, "denoised .d\n(type-2)", "#d9ead3")
    arrow(ax, 2.1, 6.9, 2.6, 6.9)
    arrow(ax, 4.7, 6.9, 5.1, 6.9)
    arrow(ax, 7.1, 6.9, 7.5, 6.9)
    ax.text(5.0, 8.0, "core: MS1 denoising (default)", ha="center", fontsize=9,
            style="italic", color="#333")
    # optional stages, branching off the survivors before write
    box(ax, 2.6, 2.6, 2.0, 1.3, "MS/MS denoise\n(dda/dia)", "#fff2cc", 8)
    box(ax, 4.9, 2.6, 2.0, 1.3, "watershed\ncentroider", "#fce5cd", 8)
    box(ax, 7.2, 2.6, 2.0, 1.3, "box-centroid\nconsolidation", "#fce5cd", 8)
    ax.text(1.0, 3.25, "optional\nstages", ha="center", va="center", fontsize=9,
            style="italic", color="#666")
    for bx in (3.6, 5.9, 8.2):
        arrow(ax, bx, 4.6, bx, 3.9)
    ax.text(6.1, 5.1, "applied to the same (MS1) frames before re-encoding",
            ha="center", fontsize=7.5, color="#666")


def main() -> int:
    fig = plt.figure(figsize=(13, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.6], wspace=0.18)
    concept(fig.add_subplot(gs[0]))
    pipeline(fig.add_subplot(gs[1]))
    fig.suptitle("dnoise: ion-mobility–aware denoising of timsTOF MS1 frames", fontsize=13)
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / "fig1_overview.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
