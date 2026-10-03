# SynTangle rule-set changelog

This file records intentional changes to the normative problem definition in `RULES.md`.

Rule changes must be explicit. A rule amendment must:

1. modify `RULES.md`;
2. update this changelog;
3. include `[RULES-CHANGE]` in the pull-request title/body or commit message;
4. state why the biological or mathematical definition of the problem is changing.

## 2026-10-02 — Initial rule set

Established R1–R16.

The initial formulation fixes the central invariant that within-chromosome genomic order is immutable. Legal layout operations are restricted to whole-chromosome movement and whole-chromosome reversal. It also establishes decomposition before optimization, propagation of forced states, separation of intrinsic from layout tangledness, and auditability of all transformations.


## 2026-10-03 — R17 monotone search-space reduction

Added R17 to make an implementation principle explicit in the normative solver definition: optimization stages must consume the reduced decision space produced by earlier stages rather than expanding previously contracted or eliminated possibilities again.

This does **not** permit biological evidence to be discarded. Full observed homology/block data remain available for scoring, reconstruction, and audit. The change constrains only the treatment of optimization variables: states proven infeasible, forced, equivalent, dominated by a valid lower bound, or unable to improve the incumbent must remain eliminated downstream.

Motivation: the Stage 12 branch-and-bound solver correctly propagated hard orientation equations into a compact free-flip basis, but then expanded that basis back into all complete orientation assignments before ordering search. R17 rules out that kind of combinatorial re-expansion and requires subsequent stages to solve only residual uncertainty.
