# Experimental ordering-only backbone

This opt-in helper generates proposals for an explicitly different ordering
objective. It does not change native evidence or redefine SynTangle's canonical
all-link crossing score. It is not a lossless kernel, a proof exclusion, or a
global optimality claim. The canonical optimiser and default inputs are unchanged.

`select_links(rows, min_support=3, min_run=3)` counts DISTINCT BUSCO IDs for each
adjacent-species chromosome pair, so duplicated copies cannot inflate support.
A retained link must also belong to consecutive same-partner runs containing
at least three distinct anchors in native order on BOTH chromosomes. Ambiguous
chromosome-partner assignments break runs. No orientation/collinearity or bp-span
claim is made; these are short-run filters by anchor count, not physical length.
Thresholds are exploratory parameters, not a previous agreed setting.

The mask audits every original link, including links not used for ordering. A
caller must preserve the full original tables, all chromosomes and coordinates,
then restore every link when plotting and calculating the full score. Isolated
chromosomes remain present. Ordering based on the mask can increase the full
score; report both objectives separately, without claiming canonical improvement.

The Pieridae experiment freezes the helper from this feature branch and uses
existing local-search whole-chromosome moves with one start, at most 50 improving
steps and conditional order DP capped at 10 chromosomes. This is a bounded-work
heuristic, not a wall-clock guarantee. It records the mask, candidate layout,
filtered solver diagnostics and full-link before/after scores. No proof
exclusion ledger from another running search is imported or modified.
