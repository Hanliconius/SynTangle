# SynTangle normative rules

> **NORMATIVE DOCUMENT**
>
> These rules define the SynTangle optimization problem. Ordinary implementation work must not modify them.
>
> Any change to this file must be an explicit rule-set amendment, must also update `RULES_CHANGELOG.md`, and must carry the marker `[RULES-CHANGE]` in the pull-request title/body or commit message. The repository guard checks this mechanically.
>
> Implementation convenience is never sufficient reason to violate a rule.

## Biological invariants

**R1 — Within-chromosome genomic order is immutable.**  
The input assembly defines the true linear order of positions, loci, anchors, and syntenic blocks along each chromosome. SynTangle must preserve that order exactly.

**R2 — Whole chromosomes are movable units.**  
A chromosome may change its position in the displayed chromosome order only as a complete unit.

**R3 — Whole chromosomes may reverse.**  
A chromosome may be displayed in either native orientation or complete reverse orientation. Reversal acts on the entire chromosome and deterministically reverses all internal coordinates/order relationships.

**R4 — Subchromosomal reordering is forbidden.**  
No optimizer, spectral method, heuristic, or visualization routine may independently move, swap, reverse, or reorder a subchromosomal block to reduce tangledness.

**R5 — Observed homology is input evidence, not an optimization variable.**  
Optimization may choose legal display states; it may not invent, delete, or reassign homology merely to improve the layout. Uncertain homology may carry explicit confidence/weight, but uncertainty must remain visible in the model.

**R6 — Missing or ambiguous data must not be coerced into false certainty.**  
Absent blocks, duplications, ambiguous mappings, assembly gaps, and low-confidence relationships must be represented explicitly rather than silently converted into one-to-one correspondence.

## Decomposition and admissible search space

**R7 — Disconnected chromosome-homology components are independent.**  
Components with no homology connection must be solved independently. Their relative global order is representationally equivalent provided a consistent component order is used across species.

**R8 — Forced states are propagated before optimization.**  
Orientation, precedence, contiguity, and component constraints that are logically determined by the data must be solved and removed from the free search space before combinatorial optimization begins.

**R9 — Only biologically admissible states may enter the optimizer.**  
The optimizer searches only whole-chromosome permutations and whole-chromosome reversals that satisfy all hard constraints derived from the input.

**R10 — True chromosome order outranks aesthetic improvement.**  
If reducing crossings would require violating extant chromosome order, that crossing is intrinsic under the SynTangle model and must remain.

**R11 — Solved/contracted structure must remain exactly recoverable.**  
Any contraction into supernodes, components, chains, or kernels must preserve enough information to reconstruct the complete legal chromosome layout without loss.

## Objective and interpretation

**R12 — Intrinsic tangledness and layout tangledness are distinct.**  
The best achievable score under all biological constraints describes residual/intrinsic tangledness. Excess tangledness of a particular drawing is measured relative to that constrained optimum/best-known solution.

**R13 — Optimization must not manufacture biological simplicity.**  
A lower visual score is valid only when obtained through legal whole-chromosome moves and flips.

**R14 — Metrics remain interpretable.**  
Structural complexity, orientation frustration, crossing tangledness, and search-space complexity are different quantities and must not be silently collapsed into one score without an explicit, documented model.

## Reproducibility

**R15 — Every transformation is auditable.**  
The final result must record the input state, every whole-chromosome move/flip required to obtain the output state, the constraints applied, and the objective values before and after optimization.

**R16 — Rule changes require explicit amendment.**  
Changing, removing, weakening, or adding a normative rule requires an intentional rules amendment recorded in `RULES_CHANGELOG.md`. Code must conform to the rules, not redefine them.
