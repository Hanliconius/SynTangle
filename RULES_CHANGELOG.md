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
