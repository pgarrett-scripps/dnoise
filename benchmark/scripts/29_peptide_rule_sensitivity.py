#!/usr/bin/env python3
"""Sensitivity of the ddaPASEF quantified-protein coverage gain to the
peptides-per-protein rule (reviewer ask: could the post-hoc two-peptide filter
interact with denoising?). Recomputes the aggregate quantified-protein count for
original vs MS1-denoised arms under the 1-, 2-, and 3-peptide rules, both gradients.

A protein is "quantified" when it clears the peptide-count rule and is present in
>=2 of the 6 replicates of each of the two compared conditions (A and B), matching
the main-text two-peptide definition.

Prints the table and is the source for SI Table S14 (tab:peptide-rule).
Usage: uv run scripts/29_peptide_rule_sensitivity.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _metrics as M  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def quant_count(arm_dir: Path, min_pep: int) -> int:
    M.MIN_PEPTIDES_PER_PROTEIN = min_pep
    prot, fc = M._rollup(arm_dir)
    if prot is None:
        return 0
    a = [c for c in fc if M.condition_of(c) == "A"]
    b = [c for c in fc if M.condition_of(c) == "B"]
    ok = (prot[a].notna().sum(axis=1) >= 2) & (prot[b].notna().sum(axis=1) >= 2)
    return int(ok.sum())


def main() -> int:
    for g in ["dda_5min", "dda_15min"]:
        print(g)
        for mp in (1, 2, 3):
            o = quant_count(ROOT / "results" / g / "original", mp)
            d = quant_count(ROOT / "results" / g / "denoised", mp)
            print(f"  {mp}-peptide rule: original {o}, MS1 {d}, "
                  f"Δ {d - o:+d} ({100 * (d - o) / o:+.1f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
