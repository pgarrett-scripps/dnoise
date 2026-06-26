# TODO — before JPR submission

## ✅ RESOLVED (2026-06-14): DIA library mismatch fixed — both gradients now FULL proteome
The two diaPASEF gradients had used DIFFERENT DIA-NN libraries (5-min = restricted
DDA-10%-FDR allowlist hybrid_diann_dda10; 15-min = full hybrid_diann), so absolute
counts were not comparable. FIX TAKEN: re-ran the 18 dia_5min runs x3 arms against
the FULL proteome library (matching 15-min). Now BOTH gradients use hybrid_diann
(31,390 isoforms / 41,867 groups / 4.84M precursors) — directly comparable.
  - Old restricted dia_5min results backed up at results/dia_5min_restricted_bak.
  - fig4_dia.png regenerated (coverage panel now full-lib for both gradients).
  - SI §S7: dropped the restricted-allowlist rationale; @tab:dia-2grad is the new
    two-gradient (7-col) DIA IDs/LFQ table; prose rewritten for both gradients.
  - Main @fig:dia caption: "restricted proteome" line replaced with "same full
    predicted proteome library".
  - paper.pdf + supplementary.pdf recompile clean; no leftover restriction text.
New dia_5min full-lib counts (orig/MS1/+MS/MS): precursors 59,343/59,592/57,137;
protein groups 9,537/9,537/9,296; quantified 8,968/8,926/8,673.

## Other pre-submission blockers (from review 2026-06-12)
- Front matter missing (JPR/ACS): co-authors, affiliation (Dept [TODO]), ORCID,
  Author Contributions, Funding, Conflict of Interest, Acknowledgements,
  Associated Content heading, designated TOC graphic.
- Citations: DIA-NN never cited (add Demichev 2020); timsrust + PXD070049 carry
  TODO notes in references.bib; intro "existing approaches" cites nothing; only 5
  refs total (thin for a full article).
- Bruker Minesweeper comparison: confirm permission/attribution (footnote flag).
- Two [TODO: exact CPU model] placeholders (paper Performance + SI S8 Hardware).
- 15-min CV framing: 0.059 -> 0.073 called "essentially unchanged" — soften.

---

overall the paper looks good.

figure 1 is great! fantastic visualization

figure2 also is good, but I want it to shwo compression for 5/15/DDA/DIA so eventaull it wuill show 4 different file/experimnt types.

figure 3 is ok.... but pretty weak. Will thik about how to make this better.

figure 4 shows identifications psm/peptide/protein/quantified. Eventaully i want this to shwo data for the 4 differetn file/experiemtns. Similar to figure 2. 

figure 5 - shows the lfq accuracy. this is good but should eventaully show for all 4 file/experiemnt types.

figure 6 - add in ms1 watershed centroider algorithm. This shoudl only be for dda 5 min. and should show the compression, IDs, and LFQ asccuracy. this one will only be run for the small testing ddatset.