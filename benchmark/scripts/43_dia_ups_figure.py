#!/usr/bin/env python3
# /// script
# dependencies = ["pandas", "pyarrow", "matplotlib"]
# ///
"""SI figure for the UPS2 / timsTOF Pro 2 second-instrument generalization arm,
showing BOTH acquisition types (ddaPASEF via Sage, diaPASEF via DIA-NN) together.

2x2 grid, both gradients in each panel (solid 30 SPD, hatched 15 SPD):
  (a) ddaPASEF: UPS proteins quantified per arm   (b) ddaPASEF: per-peptide LFQ distortion
  (c) diaPASEF: UPS proteins quantified per arm   (d) diaPASEF: per-precursor LFQ distortion

DDA arms read Sage lfq.tsv (q<=0.01) exactly as 38/39_ups_*.py; DIA arms read the
DIA-NN report.parquet (Q.Value<=0.01) exactly as 42_dia_ups_compare.py, so the
figure matches @tab:ups-ids and @tab:ups-fid by construction. Writes
paper/figures/si_ups.png.

Usage: uv run scripts/43_dia_ups_figure.py
"""

from __future__ import annotations

import csv
import math
import re
import statistics as st
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
                     "xtick.labelsize": 9.5, "ytick.labelsize": 10, "legend.fontsize": 9.5})

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT.parent / "paper" / "figures" / "si_ups.png"
FDR = 0.01
GRADS = [("30spd", "30 SPD (44 min)"), ("15spd", "15 SPD (88 min)")]
UPS_FASTA = ROOT / "data" / "fasta" / "ups_ecoli.fasta"

CONC = {r["accession"]: float(r["oncolumn_fmol"])
        for r in csv.DictReader(open(ROOT / "data/meta/ups2_concentrations.tsv"), delimiter="\t")}


def acc2sp(fasta):
    m = {}
    for line in open(fasta):
        if line.startswith(">"):
            p = line[1:].split("|")
            if len(p) >= 2:
                m[p[1]] = line[1:].split("_", 1)[0]
    return m


ACC2SP = acc2sp(UPS_FASTA)
ARMS = ["original", "denoised", "intensity", "msms", "intensity_msms"]
LABELS = ["original", "streak\nMS1", "intens.\nMS1", "streak\n+MS/MS", "intens.\n+MS/MS"]
COL = {"original": "#999999", "denoised": "#E69F00", "intensity": "#56B4E9",
       "msms": "#D55E00", "intensity_msms": "#0072B2"}
FID_PAIRS = [("denoised", "MS1: streak"), ("intensity", "MS1: intens."),
             ("msms", "+MS/MS: streak"), ("intensity_msms", "+MS/MS: intens.")]


def aligned(a, ref, keys):
    lr = [math.log2(a[k] / ref[k]) for k in keys]
    s = st.median(lr)
    return st.median([abs(x - s) for x in lr])


# ---------------- DDA (Sage lfq.tsv) ----------------

def acc_of(entry):
    m = re.search(r"sp\|(?:Cont_)?([A-Z0-9]+)\|", entry)
    return m.group(1) if m else None


def is_ups_proteins(proteins):
    return any((m := re.search(r"sp\|([A-Z0-9]+)\|", e)) and m.group(1) in CONC
               for e in proteins.split(";"))


def dda_lfq(grad, arm):
    f = ROOT / "results" / f"dda_ups_{grad}" / arm / "lfq.tsv"
    if not f.is_file():
        return None
    return list(csv.DictReader(open(f), delimiter="\t"))


def dda_ups_proteins(grad, arm):
    rows = dda_lfq(grad, arm)
    if rows is None:
        return 0
    icol = [c for c in rows[0] if c.endswith(".d")][0]
    ups_pep = defaultdict(list)
    for r in rows:
        if float(r["q_value"]) > FDR:
            continue
        entries = r["proteins"].split(";")
        accs = {acc_of(e) for e in entries}
        ups_hits = {a for a in accs if a in CONC}
        if len(ups_hits) == 1 and len(entries) == 1:
            v = r[icol]
            if v not in ("", "NaN") and float(v) > 0:
                ups_pep[next(iter(ups_hits))].append(float(v))
    return sum(1 for v in ups_pep.values() if len(v) >= 2)


def dda_pep(grad, arm):
    rows = dda_lfq(grad, arm)
    if rows is None:
        return {}
    icol = [c for c in rows[0] if c.endswith(".d")][0]
    out = {}
    for r in rows:
        if float(r["q_value"]) > FDR or not is_ups_proteins(r["proteins"]):
            continue
        v = r[icol]
        if v not in ("", "NaN") and float(v) > 0:
            out[r["peptide"]] = out.get(r["peptide"], 0.0) + float(v)
    return out


def dda_data(grad):
    ids = {a: dda_ups_proteins(grad, a) for a in ARMS}
    peps = {a: dda_pep(grad, a) for a in ARMS}
    o = peps["original"]
    k1 = set(o) & set(peps["denoised"]) & set(peps["intensity"])
    k2 = set(o) & set(peps["msms"]) & set(peps["intensity_msms"])
    keys = {"denoised": k1, "intensity": k1, "msms": k2, "intensity_msms": k2}
    fid = {a: aligned(peps[a], o, keys[a]) for a, _ in FID_PAIRS}
    return ids, fid


# ---------------- DIA (DIA-NN report.parquet) ----------------

def accs(g):
    return {a.split("-", 1)[0].strip() for a in str(g).split(";") if a.strip()}


def ups_acc(g):
    a = accs(g); h = {x for x in a if x in CONC}
    return next(iter(h)) if len(h) == 1 and len(a) == 1 else None


def dia_report(grad, arm):
    f = ROOT / "results" / f"dia_ups_{grad}" / arm / "report.parquet"
    if not f.is_file():
        return None
    df = pd.read_parquet(f)
    df = df[df["Q.Value"] <= FDR]
    if "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]
    return df


def dia_ups_proteins(grad, arm):
    df = dia_report(grad, arm)
    if df is None:
        return 0
    npep = df.groupby("Protein.Group")["Stripped.Sequence"].nunique()
    keep = set(npep[npep >= 2].index)
    return len({ups_acc(g) for g in keep if ups_acc(g) is not None})


def dia_prec(grad, arm):
    df = dia_report(grad, arm)
    if df is None:
        return {}
    out = {}
    for _, r in df.iterrows():
        if ups_acc(str(r["Protein.Group"])) is None:
            continue
        v = r["Precursor.Quantity"]
        if pd.notna(v) and float(v) > 0:
            out[str(r["Precursor.Id"])] = out.get(str(r["Precursor.Id"]), 0.0) + float(v)
    return out


def dia_data(grad):
    ids = {a: dia_ups_proteins(grad, a) for a in ARMS}
    precs = {a: dia_prec(grad, a) for a in ARMS}
    o = precs["original"]
    k1 = set(o) & set(precs["denoised"]) & set(precs["intensity"])
    k2 = set(o) & set(precs["msms"]) & set(precs["intensity_msms"])
    keys = {"denoised": k1, "intensity": k1, "msms": k2, "intensity_msms": k2}
    fid = {a: aligned(precs[a], o, keys[a]) for a, _ in FID_PAIRS}
    return ids, fid


# ---------------- plot ----------------

def panel_ids(ax, data_by_grad, title, ylab):
    x = np.arange(len(ARMS)); w = 0.38
    for j, (g, glab) in enumerate(GRADS):
        ids = data_by_grad[g][0]
        vals = [ids[a] for a in ARMS]
        bars = ax.bar(x + (j - 0.5) * w, vals, w, color=[COL[a] for a in ARMS],
                      alpha=1.0 if j == 0 else 0.6, edgecolor="black", linewidth=0.5,
                      hatch="" if j == 0 else "//")
        ax.bar_label(bars, fmt="%d", fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels(LABELS)
    ax.set_ylabel(ylab); ax.set_title(title)
    ax.set_ylim(0, max(max(d[0].values()) for d in data_by_grad.values()) * 1.18)


def panel_fid(ax, data_by_grad, title):
    xb = np.arange(len(FID_PAIRS)); w = 0.38
    for j, (g, glab) in enumerate(GRADS):
        fid = data_by_grad[g][1]
        vals = [fid[a] for a, _ in FID_PAIRS]
        bars = ax.bar(xb + (j - 0.5) * w, vals, w, color=[COL[a] for a, _ in FID_PAIRS],
                      alpha=1.0 if j == 0 else 0.6, edgecolor="black", linewidth=0.5,
                      hatch="" if j == 0 else "//")
        ax.bar_label(bars, fmt="%.3f", fontsize=7.5)
    ax.set_xticks(xb); ax.set_xticklabels([lab for _, lab in FID_PAIRS], rotation=12)
    ax.set_ylabel("median $|log_2(\\mathrm{arm}/\\mathrm{orig.})|$ (aligned)")
    ax.set_title(title)


dda = {g: dda_data(g) for g, _ in GRADS}
dia = {g: dia_data(g) for g, _ in GRADS}

fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))
panel_ids(axes[0, 0], dda, "(a) ddaPASEF: UPS proteins (1% FDR, Sage)", "UPS proteins ($\\geq$2 pep.)")
panel_fid(axes[0, 1], dda, "(b) ddaPASEF: per-peptide LFQ distortion (lower = truer)")
panel_ids(axes[1, 0], dia, "(c) diaPASEF: UPS proteins (1% prec. FDR, DIA-NN)", "UPS protein groups ($\\geq$2 pep.)")
panel_fid(axes[1, 1], dia, "(d) diaPASEF: per-precursor LFQ distortion (lower = truer)")

fig.legend(handles=[Patch(facecolor="white", edgecolor="black", label=GRADS[0][1]),
                    Patch(facecolor="white", edgecolor="black", hatch="//", label=GRADS[1][1])],
           loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0), frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(FIG, dpi=150)
print(f"wrote {FIG}")
