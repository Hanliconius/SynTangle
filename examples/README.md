# Synthetic examples / future executable tests

These examples are specifications first. Each should later become a small machine-readable fixture with expected decomposition, legal state space, and metrics.

## Case 1 — Perfect 1:1:1 correspondence

Three species, 30 chromosomes each:

```text
Sp1: 1  2  3 ... 30
Sp2: 1  2  3 ... 30
Sp3: 1  2  3 ... 30
```

Each chromosome is homologous only to the same-numbered chromosome in the other species.

Expected:

- 30 independent chromosome-homology components;
- no cross-component combinatorial optimization;
- orientation solved independently within each component;
- a consistent common component order has zero layout crossings for perfectly collinear input;
- arbitrary independent reordering of chromosome interiors is forbidden.

## Case 2 — Fusion chain without a cycle

Homology groups of interest: 1, 9, 18, 30.

```text
Species 1: one chromosome contains 1 + 30
Species 2: one chromosome contains 30 + 18
           one chromosome contains 9 + 1
Species 3: 1, 9, 18, 30 are separate
```

At the coarse incidence level this forms a connected but tree-like structure.

Expected:

- groups 1, 9, 18, 30 belong to one nontrivial component;
- incidence cycle rank = 0 for the fusion chain itself;
- leaves/bridges should remove most or all combinatorial work;
- observed order within every fused chromosome is immutable;
- only whole fused chromosomes may reverse.

## Case 3 — Close the fusion cycle

Add a species/chromosome containing:

```text
9 + 18
```

Expected:

- the previous tree gains an independent cycle;
- the 2-core/biconnected cyclic kernel is non-empty;
- optimization is restricted to that kernel rather than the rest of the genome.

## Case 4 — Balanced orientation cycle

Construct a signed cycle whose XOR sum is 0.

Expected:

- frustration = 0;
- choose one whole-chromosome orientation and propagate the rest;
- no (2^n) orientation search.

## Case 5 — Frustrated orientation cycle

Construct a signed cycle whose XOR sum is 1.

Expected:

- inconsistency detected;
- conflicting core isolated;
- no subchromosomal reversal is permitted as a “fix”.

## Case 6 — Forbidden pretty solution

A chromosome has true internal order:

```text
A — B — C — D
```

Another species creates a crossing pattern that would look cleaner if B and C were swapped.

Expected:

- legal states include `A-B-C-D` and whole-chromosome reverse `D-C-B-A`;
- `A-C-B-D` is never generated;
- residual crossings required by this constraint count as intrinsic under the model.

## Case 7 — Layout noise versus biological complexity

Start from a simple zero-crossing biological correspondence, then scramble only the displayed whole-chromosome order.

Expected:

- structural/cycle metrics unchanged;
- intrinsic/best-achievable tangledness unchanged;
- current layout tangledness increases;
- solver restores a zero-excess layout using legal whole-chromosome moves.
