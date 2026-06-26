#!/usr/bin/env python3
# /// script
# dependencies = ["pandas", "pyarrow"]
# ///
"""Streak-vs-intensity comparison on the UPS2 / timsTOF Pro 2 *diaPASEF* data
(the DIA-NN analog of 38_ups_compare.py + 39_ups_lfq_fidelity.py, which did the
Sage/DDA arms). Single-injection UPS-only files (replicates arrive later), so we
report IDs, dynamic-range linearity, E. coli false positives, and LFQ fidelity
between arms -- not the HYE A/B accuracy that 12_analyze_dia.py computes.

FDR NOTE: a UPS-only proteome (48 proteins) is far too small for DIA-NN's
protein-GROUP global q-value to be estimated -- it reports 0 protein groups at
1% global FDR and the *_matrix.tsv files come out empty. The precursor-level
search is fine (~330 precursors / ~30 proteins at 1% precursor FDR). So we read
the per-precursor report.parquet, filter Q.Value<=0.01 (precursor FDR), and do
protein inference ourselves with a two-peptide rule -- exactly as we did for Sage
on the DDA side. The E. coli library entries (absent from the UPS sample) are the
entrapment: any E. coli protein passing the same rule is an empirical false +.

Arms per dataset (built/searched by 41_dia_ups_diann.sh), each searched
independently against the same predicted UPS+E.coli library:
  original  raw                                   msms            streak MS1+MS/MS
  denoised  streak MS1                            intensity_msms  matched MS1+MS2 thr
  intensity matched per-point MS1 threshold (T1)

PART A (IDs / linearity):  per arm, UPS protein groups with >=2 distinct peptides;
  per on-column concentration group; log-log fit of PG.MaxLFQ vs fmol; E. coli
  protein groups passing the same 2-peptide rule = false positives.
PART B (LFQ fidelity):     per UPS precursor, log2(arm/reference) on the common
  precursor set across the 3 compared arms, reported BOTH raw (un-normalized:
  Precursor.Quantity, the scale/bias effect normalization would remove) and after
  total-intensity normalization over the shared UPS precursors (realistic LFQ).
  global shift, median-aligned |log2| distortion, within +/-10%, r.

Usage: uv run scripts/42_dia_ups_compare.py
"""

from __future__ import annotations

import csv
import math
import statistics as st
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FDR = 0.01
DATASETS = ["dia_ups_30spd", "dia_ups_15spd"]
ARMS = ["original", "denoised", "intensity", "msms", "intensity_msms"]
UPS_FASTA = ROOT / "data" / "fasta" / "ups_ecoli.fasta"

CONC = {}
for r in csv.DictReader(open(ROOT / "data/meta/ups2_concentrations.tsv"), delimiter="\t"):
    CONC[r["accession"]] = float(r["oncolumn_fmol"])
LEVELS = sorted(set(CONC.values()), reverse=True)


def build_acc2species(fasta: Path) -> dict[str, str]:
    acc2sp = {}
    with fasta.open() as fh:
        for line in fh:
            if line.startswith(">"):
                head = line[1:]
                sp = head.split("_", 1)[0]
                parts = head.split("|")
                if len(parts) >= 2:
                    acc2sp[parts[1]] = sp
    return acc2sp


ACC2SP = build_acc2species(UPS_FASTA)


def accs(group: str) -> set[str]:
    return {a.split("-", 1)[0].strip() for a in str(group).split(";") if a.strip()}


def is_ecoli(group: str) -> bool:
    a = accs(group)
    return bool(a) and all(ACC2SP.get(x) == "ECOLI" for x in a)


def ups_acc(group: str) -> str | None:
    """The single UPS accession this group maps to, else None (no shared/ambiguous)."""
    a = accs(group)
    hits = {x for x in a if x in CONC}
    return next(iter(hits)) if len(hits) == 1 and len(a) == 1 else None


def load_report(arm: Path) -> pd.DataFrame | None:
    """Per-precursor rows at 1% precursor FDR (Q.Value<=0.01), targets only."""
    f = arm / "report.parquet"
    if not f.is_file():
        return None
    df = pd.read_parquet(f)
    if "Q.Value" not in df.columns:
        return None
    df = df[(df["Q.Value"] <= FDR)]
    if "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]
    return df


# ---------------- PART A: IDs / linearity ----------------

def linfit(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan"), float("nan"), n
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    syy = sum((y - my) ** 2 for y in ys)
    slope = sxy / sxx if sxx else float("nan")
    r2 = (sxy ** 2) / (sxx * syy) if sxx > 0 and syy > 0 else float("nan")
    return slope, r2, n


def part_a(ds: str):
    print(f"\n################  {ds}  --  IDs / linearity  ################")
    hdr = (f"{'arm':16s} {'prec':>5s} {'UPS':>4s} " +
           " ".join(f"{l:>6g}" for l in LEVELS) + f" {'slope':>6s} {'R2':>5s} {'Eco+':>5s}")
    print(hdr)
    for arm in ARMS:
        df = load_report(ROOT / "results" / ds / arm)
        if df is None:
            print(f"{arm:16s}  (no report.parquet)")
            continue
        nprec = df["Precursor.Id"].nunique()
        # distinct peptides per protein group (two-peptide rule)
        npep = df.groupby("Protein.Group")["Stripped.Sequence"].nunique()
        keep = set(npep[npep >= 2].index)
        ups_int, eco = {}, set()
        for g in keep:
            if is_ecoli(g):
                eco.add(g)
                continue
            a = ups_acc(g)
            if a is None:
                continue
            lfq = df.loc[df["Protein.Group"] == g, "PG.MaxLFQ"].dropna()
            lfq = lfq[lfq > 0]
            if len(lfq):
                ups_int[a] = ups_int.get(a, 0.0) + float(lfq.iloc[0])
        by_lvl = {l: 0 for l in LEVELS}
        xs, ys = [], []
        for a, inten in ups_int.items():
            by_lvl[CONC[a]] += 1
            xs.append(math.log10(CONC[a])); ys.append(math.log10(inten))
        slope, r2, _ = linfit(xs, ys)
        cnts = " ".join(f"{by_lvl[l]:>6d}" for l in LEVELS)
        print(f"{arm:16s} {nprec:>5d} {len(ups_int):>4d} {cnts} {slope:>6.2f} {r2:>5.2f} {len(eco):>5d}")


# ---------------- PART B: LFQ fidelity ----------------

def load_prec(arm: Path) -> dict[str, float] | None:
    """UPS precursor (Precursor.Id) -> raw Precursor.Quantity."""
    df = load_report(arm)
    if df is None:
        return None
    out = {}
    for _, row in df.iterrows():
        if ups_acc(str(row["Protein.Group"])) is None:
            continue
        v = row["Precursor.Quantity"]
        if pd.isna(v) or float(v) <= 0:
            continue
        pid = str(row["Precursor.Id"])
        out[pid] = out.get(pid, 0.0) + float(v)
    return out


def fidelity(a: dict, ref: dict, label: str, keys, normalize: bool):
    if len(keys) < 3:
        print(f"    {label:32s}: n={len(keys)} (too few)")
        return
    ratios = {k: a[k] / ref[k] for k in keys}
    if normalize:  # total-intensity normalize a to ref over the shared set
        f = sum(ref[k] for k in keys) / sum(a[k] for k in keys)
        ratios = {k: r * f for k, r in ratios.items()}
    lr = [math.log2(r) for r in ratios.values()]
    shift = st.median(lr)
    aligned = st.median([abs(x - shift) for x in lr])
    within10 = sum(1 for r in ratios.values() if 0.909 <= r <= 1.1) / len(keys)
    la = [math.log10(a[k]) for k in keys]
    lref = [math.log10(ref[k]) for k in keys]
    n = len(keys)
    mla, mlr = sum(la) / n, sum(lref) / n
    sab = sum((x - mla) * (y - mlr) for x, y in zip(la, lref))
    saa = sum((x - mla) ** 2 for x in la)
    sbb = sum((y - mlr) ** 2 for y in lref)
    r = sab / math.sqrt(saa * sbb) if saa > 0 and sbb > 0 else float("nan")
    print(f"    {label:32s}: n={n:4d}  shift(log2)={shift:+.3f}  "
          f"|log2|aligned={aligned:.3f}  within10%={within10*100:4.0f}%  r={r:.4f}")


def part_b(ds: str):
    print(f"\n################  {ds}  --  LFQ fidelity (UPS precursors)  ################")
    arms = {a: load_prec(ROOT / "results" / ds / a) for a in ARMS}
    if not arms["original"]:
        print("  (no original arm output)")
        return
    o = arms["original"]
    for norm in (False, True):
        tag = "TOTAL-INTENSITY NORMALIZED" if norm else "RAW (un-normalized)"
        print(f"  --- {tag} ---")
        if arms["denoised"] and arms["intensity"]:
            k1 = set(o) & set(arms["denoised"]) & set(arms["intensity"])
            print(f"  MS1 level (common UPS precursors in all 3 arms: n={len(k1)}):")
            fidelity(arms["denoised"], o, "streak(denoised) vs original", k1, norm)
            fidelity(arms["intensity"], o, "intensity vs original", k1, norm)
            fidelity(arms["denoised"], arms["intensity"], "streak vs intensity", k1, norm)
        if arms["msms"] and arms["intensity_msms"]:
            k2 = set(o) & set(arms["msms"]) & set(arms["intensity_msms"])
            print(f"  MS1+MS2 level (common set: n={len(k2)}):")
            fidelity(arms["msms"], o, "streak(msms) vs original", k2, norm)
            fidelity(arms["intensity_msms"], o, "intensity_msms vs original", k2, norm)
            fidelity(arms["msms"], arms["intensity_msms"], "streak vs intensity (msms)", k2, norm)


def main() -> int:
    for ds in DATASETS:
        part_a(ds)
        part_b(ds)
    print("\nPART A: prec=precursors @1% precursor FDR; UPS=quantified UPS protein groups "
          "(>=2 peptides); per-group counts at on-column fmol; slope/R2 = log10(PG.MaxLFQ) "
          "vs log10(fmol); Eco+ = E. coli false+ (entrapment).")
    print("PART B: shift<0 = intensities biased down; |log2|aligned lower = less per-precursor "
          "distortion; within10% higher = truer; r near 1 = preserved ranking. RAW shows the "
          "scale bias normalization would remove; NORMALIZED is the realistic LFQ workflow.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
