# Method-stage utility audit

This document distinguishes machinery that currently changes the optimization
from machinery that is implemented for structural analysis, validation, or
future use.

The distinction matters: SynTangle should not retain a method merely because
it is mathematically attractive. A stage should eventually justify itself by
correctness, search-space reduction, stronger bounds, better scaling, clearer
interpretation, or audit value.

| Stage | Current status | Current role |
|---|---|---|
| Sparse chromosome--homology incidence graph | **Active solver machinery** | Defines exact connected components used by exact DP, local search, and monotone branch-and-bound. |
| Connected-component factorization | **Active solver machinery** | Splits independent optimization problems under R7. |
| GF(2) orientation propagation / OrientationBasis | **Active solver machinery** | Encodes hard orientation equations and compact residual flip degrees of freedom. |
| Species-layer dynamic programming | **Active solver machinery** | Primary exact solver when transition estimates fit configured limits. |
| One-layer subset DP | **Active solver machinery** | Used in local-search coordinate optimization and branch-and-bound conditional lower bounds/order choices. |
| Local search | **Active supporting machinery** | Produces feasible incumbents for branch-and-bound; also remains a standalone best-known solver. |
| Monotone residual branch-and-bound | **Active solver machinery** | Handles cases beyond the exact layer-DP cap and can prove optimality with bounds. |
| Structural bundle projection | **Diagnostic / experimental** | Compresses repeated relationship evidence for topology summaries. It is reported in benchmarks but does not currently restrict optimizer states. |
| Bridge/articulation/2-core/biconnected decomposition | **Diagnostic / experimental** | Computed and visualized. Stage 15 still searches complete incidence components rather than only these structural kernels. |
| Fundamental cycle basis | **Diagnostic / experimental** | Describes independent graph cycles without enumerating all simple cycles; not currently used in objective or branching decisions. |
| Derived interchromosomal ordering constraints | **Infrastructure, not current hot path** | Implemented conservatively, but current benchmark inputs do not supply genuine precedence/contiguity information that would make them restrict the main solver. |
| Spectral ordering | **Deferred** | Not used by the current solver. |
| Braid-inspired moves / braid word metrics | **Deferred** | Not used by the current solver. |
| Tree-constrained ancestral event inference | **Deferred separate problem** | Explicitly downstream of stable extant-layout optimization. |

## Redundancy criterion

A stage is not called redundant merely because one benchmark does not need it.
The appropriate test is an ablation over a benchmark family:

- remove or bypass the stage;
- hold biological simulations and tangle presentations fixed;
- compare correctness/proof status;
- compare residual state counts, lower bounds, nodes/transitions, and runtime;
- inspect whether the stage supplies unique explanatory or audit information.

The paired and scaling benchmark profiles provide controlled inputs for these
ablations.

## Current graph-theory status

Graph theory remains part of the active solver through the sparse incidence
graph and exact connected-component factorization.

The more elaborate graph-topological machinery -- structural projection,
cycle rank, bridges, articulation points, 2-core, biconnected blocks, and
cycle basis -- is currently better described as measured structural machinery
than as an active accelerator. That is intentional: the project should
demonstrate that a proposed decomposition preserves the remaining crossing
objective before using it to prune search.

A future stage can test whether residual factor-graph decomposition turns some
of this topology into a proven optimization speedup. If it does not, those
pieces should remain diagnostics rather than being presented as necessary
solver steps.
