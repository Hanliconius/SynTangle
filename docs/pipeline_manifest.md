# SynTangle pipeline manifest

This is the current intended computational sequence. It is subordinate to `RULES.md`.

## 1. Normalize

Convert input synteny/homology data to the canonical chromosome and homologous-block occurrence model. Preserve genomic coordinates, orientation, confidence, ambiguity, and provenance.

## 2. Build the multispecies hypergraph

Represent each multispecies homology group once as a hyperedge. Retain chromosome membership and immutable within-chromosome coordinates on every occurrence.

## 3. Derive chromosome-level correspondence

Collapse block evidence to chromosome-scale correspondence weights for coarse decomposition and layout reasoning. Generate pairwise projections on demand rather than globally expanding every hyperedge into all pairs.

## 4. Decompose connected chromosome-homology components

Find independent components immediately. Solve them independently. Disconnected component order is representationally equivalent when applied consistently across species.

## 5. Solve whole-chromosome orientation constraints

Construct signed orientation equations (x_i\oplus x_j=b_{ij}).

- If a component is balanced, propagate orientations directly.
- If frustration is low, isolate the conflicting core rather than searching all chromosome orientations.
- Never create a subchromosomal flip.

## 6. Derive ordering constraints

Use immutable extant chromosome interiors to derive:

- precedence constraints;
- contiguous/reversible chains;
- fusion/fission-implied grouping constraints;
- relationships whose relative order is forced.

The only reversible internal object is the **entire chromosome**.

## 7. Contract forced structure

Remove solved degrees of freedom before optimization:

- trivial 1:1 homology components;
- balanced signed regions;
- forced chains;
- leaves;
- bridge-attached structures;
- any other state whose legal placement/orientation is determined.

All contractions must be exactly reversible.

## 8. Extract the hard kernel

Within each nontrivial component, compute/use:

- bridges;
- articulation points;
- 2-core;
- biconnected components;
- cycle basis and cycle rank;
- orientation frustration;
- remaining unresolved ordering relations.

The expensive solver receives only the smallest unresolved legal kernel.

## 9. Represent the admissible ordering space

Encode the remaining legal chromosome arrangements with constrained signed permutations and, where useful, PQ/PC-tree or related consecutive/partial-order machinery.

Classify relations as:

```text
FORCED      relative order/orientation is determined
EQUIVALENT  different global placements are representation-equivalent
UNRESOLVED  multiple legal choices remain and affect tangledness
```

Only UNRESOLVED choices enter crossing optimization.

## 10. Choose the smallest adequate solver

Prefer exact methods for small kernels.

Suggested hierarchy:

```text
forced solution
  → direct enumeration / dynamic programming / branch-and-bound
  → constraint-aware local search
  → spectral/hypergraph ordering as an initializer for large kernels
```

Spectral ordering is a fallback/initializer, not permission to reorder chromosome interiors.

## 11. Minimize legal visual tangledness

Optimize whole-chromosome moves/flips only.

Primary objective:

[
C=\text{weighted synteny crossings}.
]

Also evaluate:

- bundle/block crossings (B);
- positional/diagonal displacement (D);
- local maximum crossing burden, if useful.

Use local delta scoring so a candidate move updates only affected relationships.

## 12. Re-expand solved structure

Restore contracted chromosomes/components around the solved kernels without altering the optimized legal state.

## 13. Quantify the result

Keep distinct:

- intrinsic/best-achievable crossing tangledness;
- excess layout tangledness;
- cycle rank / structural complexity;
- orientation frustration;
- unresolved search-space complexity;
- largest hard-kernel size;
- optimality status.

Do not collapse these into one magic score in version 1.

## 14. Optional diagnostics

Experimental/deferred analyses may include:

- bundle-crossing perception metrics;
- multiscale tangledness by synteny resolution;
- spectral diagnostics;
- braid-like residual descriptions.

These do not override the normative state space.

## 15. Visualize

Generate all views from the same solved state:

- multispecies chromosome/ribbon layout;
- pairwise dot-plot projections;
- chromosome correspondence matrices;
- component/core diagnostics;
- metric summaries.

## 16. Validate before general implementation

Synthetic cases must verify that:

- perfect 1:1 correspondence needs no combinatorial chromosome search;
- display scrambling increases layout tangledness without changing intrinsic structure;
- whole-chromosome reversal can be corrected;
- forbidden subchromosomal reordering is never used;
- simple fusion/fission trees decompose without unnecessary optimization;
- closing a cycle creates the expected hard kernel;
- frustrated orientation cycles are detected;
- repeated starts return identical/equivalent optima where expected.
