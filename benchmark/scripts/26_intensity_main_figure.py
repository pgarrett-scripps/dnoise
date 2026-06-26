#!/usr/bin/env python3
"""Main-text figure: streak filter vs. a matched strict intensity threshold.

Promotes the intensity-baseline control (SI Section S6) to a main-text figure
(paper/figures/fig_intensity.png), showing for all three arms (original, streak
filter, intensity threshold), at both gradients:
  (A) quantified proteins, (B) quantified peptides, (C) LFQ accuracy.

Also runs a paired permutation test (sign-flip, pure numpy) on the per-run
quantified-protein counts (streak vs. intensity, 18 runs per gradient) and prints
the observed per-run difference and p-value for citation in the text.

Quant uses the shared two-peptide rule in _metrics.py. Run:
  uv run scripts/26_intensity_main_figure.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from _metrics import SPECIES, lfq_metrics, lfq_table, pair_accuracy

plt.rcParams.update({
    "font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 11,
})

ROOT = Path(__file__).resolve().parents[1]
FIGS = ROOT.parent / "paper" / "figures"
# arm label -> results subdir
ARM = {"original": "original", "streak filter": "denoised", "intensity threshold": "intensity"}
ARM_COLOR = {"original": "#0072B2", "streak filter": "#E69F00", "intensity threshold": "#D55E00"}
GRADS = [("5 min", "dda_5min"), ("15 min", "dda_15min")]
MARK = {"HUMAN": "o", "YEAST": "s", "ECOLI": "^"}


def perrun_protein_counts(arm_dir: Path) -> np.ndarray:
    import _metrics as M
    prot, fcols = M._rollup(arm_dir)
    return prot[fcols].notna().sum(axis=0).values


def paired_perm_p(a, b, iters: int = 100000, seed: int = 0) -> tuple[float, float]:
    """Two-sided paired sign-flip permutation test on the mean difference."""
    d = np.asarray(a, float) - np.asarray(b, float)
    obs = float(d.mean())
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(iters, len(d)))
    null = (signs * np.abs(d)).mean(axis=1)
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (iters + 1)
    return obs, float(p)


def main() -> int:
    # ---- collect metrics ----
    metrics = {}   # (grad_label, arm_label) -> lfq_metrics dict
    acc = {}       # (grad_label, arm_label) -> accuracy list
    for glab, ds in GRADS:
        for alab, rsub in ARM.items():
            d = ROOT / "results" / ds / rsub
            prot = lfq_table(d)
            metrics[(glab, alab)] = lfq_metrics(prot, d)
            acc[(glab, alab)] = pair_accuracy(alab, d)

    # ---- paired permutation test: streak vs intensity per-run protein counts ----
    print("Paired permutation test (per-run quantified proteins, streak vs intensity):")
    for glab, ds in GRADS:
        s = perrun_protein_counts(ROOT / "results" / ds / "denoised")
        i = perrun_protein_counts(ROOT / "results" / ds / "intensity")
        obs, p = paired_perm_p(s, i)
        print(f"  {glab}: streak mean {s.mean():.1f}, intensity mean {i.mean():.1f}, "
              f"Δ={obs:+.1f}/run, n={len(s)} runs, p={p:.2e}")

    # ---- figure: A proteins, B peptides, C accuracy ----
    arms = list(ARM)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.3))
    x = np.arange(len(GRADS))
    w = 0.8 / len(arms)
    for j, alab in enumerate(arms):
        off = (-(len(arms) - 1) / 2 + j) * w
        prot_vals = [metrics[(g, alab)]["n_quantified"] for g, _ in GRADS]
        pep_vals = [metrics[(g, alab)]["n_quant_peptide"] for g, _ in GRADS]
        b0 = ax[0].bar(x + off, prot_vals, w, label=alab, color=ARM_COLOR[alab])
        b1 = ax[1].bar(x + off, pep_vals, w, label=alab, color=ARM_COLOR[alab])
        ax[0].bar_label(b0, fmt="%d", fontsize=8, padding=1)
        ax[1].bar_label(b1, fmt="%d", fontsize=8, padding=1)
    for a, ttl, ylab in [(ax[0], "Proteins quantified", "proteins (1% FDR, ≥2 pep.)"),
                         (ax[1], "Peptides quantified", "quantified peptides")]:
        a.set_title(ttl); a.set_ylabel(ylab)
        a.set_xticks(x); a.set_xticklabels([g for g, _ in GRADS])
    ax[0].legend(fontsize=9)

    # Panel C: accuracy (observed vs expected), pooled over both gradients
    lim = (-3.5, 4.0)
    ax[2].plot(lim, lim, color="gray", ls="--", lw=1, zorder=0)
    for alab in arms:
        for glab, _ in GRADS:
            for a in acc[(glab, alab)]:
                ax[2].scatter(a["expected"], a["observed"], color=ARM_COLOR[alab],
                              marker=MARK[a["species"]], s=45, alpha=0.75)
        ax[2].scatter([], [], color=ARM_COLOR[alab], label=alab)
    ax[2].set_xlim(*lim); ax[2].set_ylim(*lim)
    ax[2].set_xlabel("expected log2 ratio"); ax[2].set_ylabel("observed median log2")
    ax[2].set_title("LFQ accuracy (○ human △ ecoli □ yeast)")
    ax[2].legend(fontsize=9)

    fig.suptitle("Streak filter vs. a matched strict intensity threshold "
                 "(equal MS1-point removal; MS/MS untouched)", fontsize=13)
    fig.tight_layout()
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / "fig_intensity.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
