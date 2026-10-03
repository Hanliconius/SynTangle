# SynTangle pipeline manifest

This document is the concise computational narrative for the current SynTangle
solver. It is subordinate to `RULES.md`.

The key distinction is between the **biological evidence representation** and
the **residual decision problem**. Biological evidence is never simplified
away merely to make the optimizer easier.

## 1. Normalize immutable biological evidence

Convert input synteny/homology data into chromosomes plus ordered homologous
block occurrences. Preserve genomic coordinates, strand/orientation evidence,
ambiguity, confidence, and provenance.

Within-chromosome order is immutable throughout the pipeline.

## 2. Represent multispecies homology as a hypergraph

A multispecies homology group is one n-ary relationship, not a collection of
independent pairwise edges.

Computationally, SynTangle stores this hypergraph through its sparse bipartite
incidence multigraph:

```text
chromosome vertices  <->  homology-group vertices
```

Each block occurrence remains a distinct incidence edge. Pairwise projections
are derived only when a solver/objective term needs them.

## 3. Factor exact incidence components

Find disconnected chromosome--homology incidence components immediately.
Under the current adjacent-layer crossing objective these components are exact
additive optimization factors.

They may therefore be solved independently. Monotone branch-and-bound can
optionally dispatch independent incidence components to separate worker
processes and recombine their proven/bounded results deterministically.

## 4. Propagate legal orientation structure

Encode hard whole-chromosome orientation relations as GF(2) equations:

```text
x_i XOR x_j = b_ij
```

A consistent component is represented compactly as one base assignment plus
free flip groups. The solver does not pre-expand all 2^k orientation states.

Contradictory hard orientation evidence remains an explicit infeasibility and
is not silently coerced.

## 5. Canonicalize equivalent component placement

Disconnected components can be placed in a common canonical order across
species without changing the biological evidence. This removes representational
noise before expensive search.

Canonicalization is an equivalence reduction, not a biological edit.

## 6. Build the residual decision/factor graph

After exact component and orientation reductions, distinguish:

- unresolved whole-component chromosome-order variables;
- unresolved GF(2) free-orientation variables;
- adjacent-species crossing factors that depend on those variables.

This residual variable/factor graph is now constructed explicitly.

If the residual factor graph disconnects, those pieces are mathematically
independent under the current objective. Stage 18 records those pieces,
objective-neutral variables, articulation variables, and a greedy min-fill
treewidth upper bound.

Current solver status:

- exact incidence components: **actively parallelizable**;
- finer disconnected residual-factor pieces: **represented and measured**;
- separator/tree-decomposition solving inside one connected residual component:
  **next solver integration step**.

No structural split is allowed to prune search unless objective independence is
proved.

## 7. Establish a feasible incumbent

Use constraint-aware local search over legal whole-chromosome operations to
obtain a good feasible layout.

Local moves may:

- reorder complete chromosomes within an incidence component;
- flip complete GF(2) free groups.

They may never edit chromosome interiors.

The incumbent supplies the upper bound used by exact/bounded search.

## 8. Choose the smallest adequate exact solver

The current hierarchy is:

```text
exact species-layer dynamic programming
        |
        | if transition/permutation space is too large
        v
monotone component branch-and-bound
        |
        | if node budget is exhausted
        v
explicit bounded best-known result
```

For a fixed orientation assignment, the crossing objective factors across
adjacent species layers, enabling dynamic programming rather than a full
Cartesian product of all species permutations.

One-layer subset DP is reused inside local search and branch-and-bound for
conditional order optimization and lower bounds.

## 9. Monotonically reduce residual uncertainty

Branch-and-bound operates only on decisions that remain unresolved.

At every stage it may:

- force a variable when one alternative cannot beat the incumbent;
- prune a branch whose valid lower bound cannot improve the incumbent;
- memoize repeated bound/order subproblems;
- terminate a component immediately when lower and upper bounds meet.

R17 requires that a removed degree of freedom never re-enter a later search
stage.

The biological evidence remains intact even when search variables disappear.

## 10. Recombine exact factors

Solutions for independent incidence components are combined in canonical
component order. Component lower bounds and upper bounds add.

A global result is reported as:

- **proven optimum** only when global lower and upper bounds agree;
- **bounded best known** otherwise, with the explicit remaining gap.

Future residual-factor and separator parallelism must preserve the same proof
accounting.

## 11. Quantify structure, search difficulty, and result separately

Keep distinct:

- initial/display crossing count;
- constrained optimum or best-known upper bound;
- lower bound and optimality gap;
- excess layout tangledness;
- raw incidence cycle rank;
- structural-bundle cycle rank;
- hard-core diagnostics;
- free orientation groups / implicit orientation states;
- residual factor-component sizes;
- residual articulation variables;
- residual min-fill treewidth upper bound;
- states/nodes actually evaluated;
- runtime.

Do not collapse these into a single opaque complexity score.

## 12. Visualize from the same solved state

Generate:

- hidden simulator-native baseline when available for validation only;
- deliberately tangled/public input;
- optimized layout;
- raw incidence graph;
- structural projection/core diagnostics;
- solver and proof diagnostics.

Visual inspection is part of validation, not decoration.

## 13. Validate by controlled forward simulation

The minimum benchmark logic is:

```text
known ancestor
   -> logged fusion/fission/inversion history
   -> extant genomes
   -> hidden native display
   -> independent whole-chromosome tangle
   -> PUBLIC benchmark input
   -> SynTangle solve
   -> reveal hidden validation information
```

The solver never receives the hidden evolutionary or tangle logs.

Use:

- tiny exact fixtures for invariant/correctness tests;
- paired mild/strong/random presentations of identical biology for presentation
  invariance;
- increasing species/chromosome/event complexity to expose scaling walls;
- visual native/tangled/optimized reports to ground-truth numerical claims.

## 14. Benchmark engineering requirements

Long benchmark runs must:

- print START/DONE progress for every case;
- checkpoint completed cases immediately;
- support resuming from the checkpoint;
- preserve deterministic seeds and case order;
- avoid increasing search caps merely to hide a scaling failure.

A bounded result at high complexity is useful evidence about the next solver
bottleneck.

## 15. Next exact decomposition layer

When the residual graph remains connected, test whether small separators make
it conditionally separable.

The intended hierarchy is:

```text
independent benchmark cases
  -> independent incidence components
  -> disconnected residual factor components
  -> small-separator / tree-decomposition subproblems
  -> parallel branch-and-bound frontier
  -> bounded/heuristic fallback only if still necessary
```

Separator or treewidth machinery becomes active only after demonstrating exact
objective factorization and preserving global lower/upper-bound proofs.

## 16. Deferred biological inference

Reconstructing ancestral fusion/fission/inversion histories is a separate,
later problem. When implemented, it will be constrained by the supplied species
tree and benchmarked against hidden simulator histories.

It must not be conflated with the primary extant-layout optimization problem.
