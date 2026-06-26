#!/usr/bin/env python3
"""SI table: per-species, per-condition-pair DIFFERENCE in median log2 ratio
between the original (undenoised) and MS1-denoised arms, with a paired bootstrap
95% CI on the difference.

Addresses the reviewer request to quantify, with uncertainty, how much MS1
denoising shifts relative quantification away from the original arm. For each
gradient x condition-pair x species we restrict to proteins quantified in BOTH
arms (the two-peptide rule, applied identically to each arm), compute the median
log2 ratio in each arm, and report the difference (MS1 - original) with a paired
percentile bootstrap CI (resampling the shared protein set). A CI straddling zero
means denoising introduces no detectable accuracy shift for that cell.

Writes paper/si/accuracy_diff.typ (a #table included by supplementary.typ).
Usage: uv run scripts/27_accuracy_diff.py [dda_5min dda_15min]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _metrics import PAIRS, SPECIES, condition_of, _rollup

ROOT = Path(__file__).resolve().parents[1]
SI = ROOT.parent / "paper" / "si"
GRAD_LABEL = {"dda_5min": "5 min", "dda_15min": "15 min"}
SP_LABEL = {"HUMAN": "human", "YEAST": "yeast", "ECOLI": "E. coli"}
N_BOOT = 2000
SEED = 0


def per_protein_ratio(prot: pd.DataFrame, fcols: list[str], a: str, b: str) -> pd.DataFrame:
    """log2(mean_a / mean_b) per (protein, species), each condition needing >=2 reps."""
    a_cols = [c for c in fcols if condition_of(c) == a]
    b_cols = [c for c in fcols if condition_of(c) == b]
    a_ok = prot[a_cols].notna().sum(axis=1) >= 2
    b_ok = prot[b_cols].notna().sum(axis=1) >= 2
    a_mean = prot[a_cols].mean(axis=1).where(a_ok)
    b_mean = prot[b_cols].mean(axis=1).where(b_ok)
    ratio = np.log2(a_mean / b_mean)
    out = pd.DataFrame({"protein": prot["protein"], "species": prot["species"],
                        "ratio": ratio})
    return out.dropna(subset=["ratio"])


def paired_diff_ci(orig: np.ndarray, ms1: np.ndarray):
    """Observed median(ms1)-median(orig) and a paired bootstrap 95% CI (rows aligned)."""
    n = len(orig)
    obs = float(np.median(ms1) - np.median(orig))
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, n, size=(N_BOOT, n))
    diffs = np.median(ms1[idx], axis=1) - np.median(orig[idx], axis=1)
    return obs, float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def rows_for(dataset: str) -> list[dict]:
    orig_prot, ofc = _rollup(ROOT / "results" / dataset / "original")
    ms1_prot, mfc = _rollup(ROOT / "results" / dataset / "denoised")
    if orig_prot is None or ms1_prot is None:
        print(f"  {dataset}: missing original/denoised rollup -- skipped", file=sys.stderr)
        return []
    out = []
    for a, b in PAIRS:
        o = per_protein_ratio(orig_prot, ofc, a, b)
        m = per_protein_ratio(ms1_prot, mfc, a, b)
        for sp in SPECIES:
            os_ = o[o["species"] == sp][["protein", "ratio"]].rename(columns={"ratio": "o"})
            ms_ = m[m["species"] == sp][["protein", "ratio"]].rename(columns={"ratio": "m"})
            shared = os_.merge(ms_, on="protein", how="inner")
            if len(shared) < 3:
                continue
            obs, lo, hi = paired_diff_ci(shared["o"].to_numpy(), shared["m"].to_numpy())
            out.append({
                "grad": GRAD_LABEL[dataset], "pair": f"{a}/{b}", "species": SP_LABEL[sp],
                "n": len(shared),
                "med_o": float(shared["o"].median()), "med_m": float(shared["m"].median()),
                "delta": obs, "lo": lo, "hi": hi,
                "straddles0": lo <= 0 <= hi,
            })
            print(f"  {dataset} {a}/{b} {sp:5s}: n={len(shared):4d} "
                  f"orig={shared['o'].median():+.3f} MS1={shared['m'].median():+.3f} "
                  f"Δ={obs:+.3f} [{lo:+.3f}, {hi:+.3f}]"
                  f"{'' if lo <= 0 <= hi else '  *'}")
    return out


def emit_typst(rows: list[dict]) -> None:
    SI.mkdir(parents=True, exist_ok=True)
    lines = [
        "#table(",
        "  columns: 8,",
        "  align: (left, center, center, center, center, center, center, center),",
        "  table.header(",
        "    [Gradient], [Pair], [Species], [$n$], [Median orig.], [Median MS1],",
        "    [$Delta$ (MS1 #sym.minus orig.)], [95% CI of $Delta$],",
        "  ),",
    ]
    for r in rows:
        pair = r["pair"].replace("/", "\\/")
        ci = f"[{r['lo']:+.3f}, {r['hi']:+.3f}]".replace("-", "#sym.minus ").replace("+", "+")
        delta = f"{r['delta']:+.3f}".replace("-", "#sym.minus ")
        mo = f"{r['med_o']:+.3f}".replace("-", "#sym.minus ")
        mm = f"{r['med_m']:+.3f}".replace("-", "#sym.minus ")
        lines.append(
            f"  [{r['grad']}], [{pair}], [_{r['species']}_], [{r['n']}], "
            f"[{mo}], [{mm}], [{delta}], [{ci}],"
        )
    lines.append(")")
    out = SI / "accuracy_diff.typ"
    out.write_text("\n".join(lines) + "\n")
    n_straddle = sum(r["straddles0"] for r in rows)
    print(f"\nwrote {out} ({len(rows)} rows; {n_straddle}/{len(rows)} CIs straddle 0)")


def main() -> int:
    datasets = sys.argv[1:] or ["dda_5min", "dda_15min"]
    rows: list[dict] = []
    for ds in datasets:
        print(f"=== {ds} ===")
        rows += rows_for(ds)
    if rows:
        emit_typst(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
