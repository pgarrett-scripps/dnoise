#!/usr/bin/env python3
"""Select the timsTOF files for one benchmark DATASET from PXD070049.

PXD070049 ("LFQ Benchmark Dataset - Generation Beta") is an HYE 3-species mix
(Human / S. cerevisiae / E. coli) acquired on many instruments. We use the
timsTOF Ultra 2 runs (type-2 .d, dnoise-readable). The deposit offers each
acquisition mode at exactly two gradients, 5 min and 15 min (nothing longer).

A "dataset" here = one acquisition group (instrument/mode/gradient/load) over a
fixed condition x replicate design. We use the mixed conditions A/B/C in six
replicates each -> 18 files per dataset.

Mixing ratios (from the dataset SDRF), w/w Human/Yeast/Ecoli:
  Condition A = 0.65 / 0.30 / 0.05
  Condition B = 0.65 / 0.15 / 0.20
  Condition C = 0.65 / 0.03 / 0.32
=> expected log2(A/B): HUMAN 0, YEAST +1, ECOLI -2 (see scripts/06_analyze.py).

Usage:  uv run scripts/01_list_files.py <dataset>
where <dataset> is one of the keys in DATASETS below (default dda_5min).
Writes config/files.<dataset>.tsv and halts if the file count is off.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.request
from pathlib import Path

PROJECT = "PXD070049"
API = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{PROJECT}/files"

CONDITIONS = ("A", "B", "C")
REPLICATES = (1, 2, 3, 4, 5, 6)

# dataset name -> PRIDE acquisition group prefix (file = <group>_Condition_<C>_REP<r>.d.zip)
DATASETS = {
    "dda_5min": "LFQ_Ultra2_PASEF_5min_50ng",
    "dda_15min": "LFQ_Ultra2_PASEF_15min_50ng",
    "dia_5min": "LFQ_Ultra2_diaPASEF_5min_50ng",
    "dia_15min": "LFQ_Ultra2_diaPASEF_15min_50ng",
}

CONFIG = Path(__file__).resolve().parents[1] / "config"


def fetch(url: str, tries: int = 8) -> list[dict]:
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=45) as r:
                return json.load(r)
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2)
    return []


def fetch_all() -> list[dict]:
    files: list[dict] = []
    page = 0
    while True:
        chunk = fetch(f"{API}?pageSize=100&page={page}")
        if not chunk:
            break
        files.extend(chunk)
        page += 1
        if page > 60:
            break
    return files


def ftp_url(f: dict) -> str:
    return next(loc["value"] for loc in f["publicFileLocations"] if "FTP" in loc["name"])


def main() -> int:
    dataset = sys.argv[1] if len(sys.argv) > 1 else "dda_5min"
    if dataset not in DATASETS:
        print(f"unknown dataset '{dataset}'. known: {', '.join(DATASETS)}", file=sys.stderr)
        return 2
    group = DATASETS[dataset]
    out = CONFIG / f"files.{dataset}.tsv"

    files = fetch_all()
    by_name = {f["fileName"]: f for f in files}

    selected = []
    for cond in CONDITIONS:
        for rep in REPLICATES:
            name = f"{group}_Condition_{cond}_REP{rep}.d.zip"
            f = by_name.get(name)
            if f is None:
                print(f"ERROR: expected file not found in {PROJECT}: {name}", file=sys.stderr)
                return 1
            selected.append((cond, f"{rep:02d}", name, f["fileSizeBytes"], ftp_url(f)))

    selected.sort()
    expected_per_cond = len(REPLICATES)
    counts = {c: sum(1 for x in selected if x[0] == c) for c in CONDITIONS}
    total = sum(s for *_, s, _ in selected)
    print(f"[{dataset}] selected {len(selected)} files from {PROJECT} "
          f"({', '.join(f'{c}:{n}' for c, n in counts.items())})")
    print(f"Total (zipped): {total/1e9:.1f} GB")

    if any(n != expected_per_cond for n in counts.values()):
        print(f"\nERROR: expected {expected_per_cond} files per condition.", file=sys.stderr)
        return 1

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("condition\treplicate\tfilename\tsize_bytes\tftp_url\n")
        for cond, rep, name, size, url in selected:
            fh.write(f"{cond}\t{rep}\t{name}\t{size}\t{url}\n")
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
