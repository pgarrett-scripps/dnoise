#!/usr/bin/env python3
"""LFQ intensity fidelity on the UPS2 / Pro 2 DDA data: do the per-peptide and
per-protein LFQ intensities change between arms? Ideally denoising leaves the
quantified intensities of real (UPS) peptides unaffected.

For each arm we read lfq.tsv (q_value<=0.01), restrict to UPS peptides, key by
peptide sequence, and compare intensities to a reference arm on the SHARED peptides:
  global shift       = median log2(arm/ref)        (0 = same scale; <0 = biased down)
  per-peptide distort= median |log2(arm/ref) - shift|   (median-aligned; lower=truer)
  within +/-10%      = fraction with 0.909 <= arm/ref <= 1.1 (raw, not aligned)
  Pearson r          = correlation of log10 intensities

Comparisons (both gradients):
  MS1 level     : denoised vs original, intensity vs original, denoised vs intensity
  MS1+MS2 level : msms vs original, intensity_msms vs original, msms vs intensity_msms

Usage: uv run scripts/39_ups_lfq_fidelity.py
"""

from __future__ import annotations

import csv
import math
import re
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FDR = 0.01
DATASETS = ["dda_ups_30spd", "dda_ups_15spd"]
UPS = {r["accession"] for r in csv.DictReader(
    open(ROOT / "data/meta/ups2_concentrations.tsv"), delimiter="\t")}


def is_ups(proteins: str) -> bool:
    return any(re.search(r"sp\|([A-Z0-9]+)\|", e) and
               re.search(r"sp\|([A-Z0-9]+)\|", e).group(1) in UPS
               for e in proteins.split(";"))


def load_pep(arm_dir: Path):
    """UPS peptide -> intensity (q<=FDR, intensity>0)."""
    f = arm_dir / "lfq.tsv"
    if not f.is_file():
        return None
    rows = list(csv.DictReader(open(f), delimiter="\t"))
    icol = [c for c in rows[0] if c.endswith(".d")][0]
    out = {}
    for r in rows:
        if float(r["q_value"]) > FDR or not is_ups(r["proteins"]):
            continue
        v = r[icol]
        if v in ("", "NaN"):
            continue
        v = float(v)
        if v > 0:
            out[r["peptide"]] = out.get(r["peptide"], 0.0) + v
    return out


def compare(a: dict, ref: dict, label: str, keys=None):
    keys = (set(a) & set(ref)) if keys is None else keys
    if len(keys) < 3:
        print(f"    {label:28s}: n={len(keys)} (too few)")
        return
    lr = [math.log2(a[k] / ref[k]) for k in keys]
    shift = st.median(lr)
    aligned = [abs(x - shift) for x in lr]
    within10 = sum(1 for k in keys if 0.909 <= a[k] / ref[k] <= 1.1) / len(keys)
    la = [math.log10(a[k]) for k in keys]
    lref = [math.log10(ref[k]) for k in keys]
    n = len(keys)
    mla, mlr = sum(la) / n, sum(lref) / n
    sab = sum((x - mla) * (y - mlr) for x, y in zip(la, lref))
    saa = sum((x - mla) ** 2 for x in la)
    sbb = sum((y - mlr) ** 2 for y in lref)
    r = sab / math.sqrt(saa * sbb) if saa > 0 and sbb > 0 else float("nan")
    print(f"    {label:28s}: n={n:4d}  shift(log2)={shift:+.3f}  "
          f"|log2|aligned={st.median(aligned):.3f}  within10%={within10*100:4.0f}%  r={r:.4f}")


def main() -> int:
    for ds in DATASETS:
        print(f"\n################  {ds}  ################")
        arms = {a: load_pep(ROOT / "results" / ds / a)
                for a in ["original", "denoised", "intensity", "msms", "intensity_msms"]}
        o = arms["original"]
        # common peptide set so the three comparisons reconcile exactly
        k1 = set(o) & set(arms["denoised"]) & set(arms["intensity"])
        print(f"  MS1 level (common UPS peptides quantified in all 3 arms: n={len(k1)}):")
        compare(arms["denoised"], o, "streak(denoised) vs original", k1)
        compare(arms["intensity"], o, "intensity vs original", k1)
        compare(arms["denoised"], arms["intensity"], "streak vs intensity", k1)
        k2 = set(o) & set(arms["msms"]) & set(arms["intensity_msms"])
        print(f"  MS1+MS2 level (common set: n={len(k2)}):")
        compare(arms["msms"], o, "streak(msms) vs original", k2)
        compare(arms["intensity_msms"], o, "intensity_msms vs original", k2)
        compare(arms["msms"], arms["intensity_msms"], "streak vs intensity (msms)", k2)
    print("\nshift<0 = intensities biased downward; |log2|aligned lower = less per-peptide "
          "distortion; within10% higher = truer; r near 1 = preserved.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
