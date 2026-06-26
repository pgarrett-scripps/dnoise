#!/usr/bin/env python3
"""Streak-vs-intensity comparison on the UPS2 / timsTOF Pro 2 DDA data, at the MS1
and MS1+MS2 levels, at matched data reduction.

Arms per dataset:
  original         raw
  denoised         streak MS1 (vertical+halo+polygon)        | MS1 level: vs intensity
  intensity        matched per-point MS1 threshold (T1)      |
  msms             streak MS1 + streak MS/MS                 | MS1+MS2 level: vs intensity_msms
  intensity_msms   matched per-point MS1 (T1) + MS2 (T2)     |

For each arm: read lfq.tsv (Sage), keep q_value<=0.01, assign peptides uniquely to
a UPS protein (no shared/ambiguous), require >=2 peptides (two-peptide rule), and
sum peptide intensities -> protein intensity. Then per UPS concentration group
(on-column fmol from data/meta/ups2_concentrations.tsv) report detected UPS
proteins and a log-log linearity fit (slope, R^2). E. coli proteins passing the
same rule are false positives (FDR sanity).

Usage: uv run scripts/38_ups_compare.py
"""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FDR = 0.01
DATASETS = ["dda_ups_30spd", "dda_ups_15spd"]
ARMS = ["original", "denoised", "intensity", "msms", "intensity_msms"]

# accession -> on-column fmol
CONC = {}
for r in csv.DictReader(open(ROOT / "data/meta/ups2_concentrations.tsv"), delimiter="\t"):
    CONC[r["accession"]] = float(r["oncolumn_fmol"])
LEVELS = sorted(set(CONC.values()), reverse=True)


def acc_of(entry: str):
    m = re.search(r"sp\|(?:Cont_)?([A-Z0-9]+)\|", entry)
    return m.group(1) if m else None


def parse_arm(arm_dir: Path):
    """Return {ups_accession: protein_intensity} (>=2 unique peptides) and E. coli
    false-positive protein count, from lfq.tsv at 1% q."""
    f = arm_dir / "lfq.tsv"
    if not f.is_file():
        return None, None
    rows = list(csv.DictReader(open(f), delimiter="\t"))
    icol = [c for c in rows[0] if c.endswith(".d")][0]
    ups_pep = defaultdict(list)   # acc -> [intensity per unique peptide]
    eco = set()
    for r in rows:
        if float(r["q_value"]) > FDR:
            continue
        entries = r["proteins"].split(";")
        accs = {acc_of(e) for e in entries}
        ups_hits = {a for a in accs if a in CONC}
        if len(entries) == 1 and any("ECOLI_" in e for e in entries):
            eco.add(accs.pop())
            continue
        # assign to a UPS protein only if it maps to exactly one UPS accession
        if len(ups_hits) == 1 and len(entries) == 1:
            inten = float(r[icol]) if r[icol] not in ("", "NaN") else 0.0
            if inten > 0:
                ups_pep[next(iter(ups_hits))].append(inten)
    prot = {a: sum(v) for a, v in ups_pep.items() if len(v) >= 2}
    return prot, len(eco)


def linfit(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan"), float("nan"), n
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    syy = sum((y - my) ** 2 for y in ys)
    slope = sxy / sxx
    r2 = (sxy ** 2) / (sxx * syy) if sxx > 0 and syy > 0 else float("nan")
    return slope, r2, n


def main() -> int:
    for ds in DATASETS:
        print(f"\n################  {ds}  ################")
        hdr = f"{'arm':16s} {'UPS':>4s} " + " ".join(f"{l:>6g}" for l in LEVELS) + f" {'slope':>6s} {'R2':>5s} {'Eco+':>5s}"
        print(hdr)
        arms_prot = {}
        for arm in ARMS:
            prot, eco = parse_arm(ROOT / "results" / ds / arm)
            if prot is None:
                print(f"{arm:16s}  (no lfq.tsv)")
                continue
            arms_prot[arm] = prot
            by_lvl = {l: 0 for l in LEVELS}
            xs, ys = [], []
            for a, inten in prot.items():
                lvl = CONC[a]
                by_lvl[lvl] += 1
                xs.append(math.log10(lvl)); ys.append(math.log10(inten))
            slope, r2, n = linfit(xs, ys)
            cnts = " ".join(f"{by_lvl[l]:>6d}" for l in LEVELS)
            print(f"{arm:16s} {len(prot):>4d} {cnts} {slope:>6.2f} {r2:>5.2f} {eco:>5d}")
        # explicit head-to-head deltas
        def npro(a): return len(arms_prot.get(a, {}))
        if "denoised" in arms_prot and "intensity" in arms_prot:
            print(f"  MS1 level   : streak={npro('denoised')}  intensity={npro('intensity')}  "
                  f"(IDs identical by construction; compare linearity/retention above)")
        if "msms" in arms_prot and "intensity_msms" in arms_prot:
            print(f"  MS1+MS2 lvl : streak={npro('msms')}  intensity={npro('intensity_msms')}  "
                  f"Δ={npro('msms')-npro('intensity_msms'):+d} UPS proteins (matched MS2 removal)")
    print("\ncolumns: UPS=quantified UPS proteins (>=2 peptides); per-group counts at "
          "on-column fmol; slope/R2 = log10(intensity) vs log10(fmol) linearity; "
          "Eco+ = E. coli false-positive proteins.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
