# Method-stage utility audit

This document distinguishes machinery that currently changes optimization from
machinery that is diagnostic, preparatory, or deferred.

A method belongs in the claimed core only if it contributes to correctness,
search-space reduction, stronger bounds, better scaling, interpretation, or
auditability. Mathematical appeal alone is not sufficient.

| Stage | Current status | Current role |
|---|---|---|
| Multispecies homology hypergraph / incidence multigraph | **Active solver representation** | Preserves n-ary homology through an exact sparse chromosome↔homology incidence representation. |
| Exact incidence connected components | **Active solver machinery** | Factor the objective into independent problems under R7. Monotone branch-and-bound can solve these components in separate worker processes. |
| Component canonicalization | **Active solver machinery** | Removes representation-equivalent global component placement before expensive search. |
| GF(2) orientation propagation / OrientationBasis | **Active solver machinery** | Encodes hard orientation equations and compact free whole-chromosome flip degrees of freedom. |
| Exact species-layer dynamic programming | **Active solver machinery** | Primary exact solver when transition/permutation estimates fit configured limits. Currently serial across incidence components. |
| One-layer subset DP | **Active solver machinery** | Used by local search and branch-and-bound for conditional chromosome-order optimization and bounds. |
| Constraint-aware local search | **Active supporting machinery** | Produces legal feasible incumbents and therefore useful upper bounds. |
| Monotone component branch-and-bound | **Active solver machinery** | Handles cases beyond exact layer-DP limits; reports explicit LB/UB/gap and can parallelize independent incidence components. |
| Residual variable/factor graph | **Active diagnostic representation; next solver factorization layer** | Explicitly records unresolved order/orientation variables and the adjacent-layer crossing factors that couple them. Disconnected residual pieces, neutral variables, articulation variables, and min-fill treewidth upper bounds are now measured. |
| Residual disconnected-piece solving | **Next integration step** | Exact factorization is conceptually available when the residual factor graph disconnects, but the current generic solver has not yet been rewritten to dispatch those finer pieces independently. |
| Small separators / tree decomposition | **Experimental next step** | Candidate mechanism for conditionally splitting a single connected residual problem; may become active only after proof-preserving recombination is implemented. |
| Structural bundle projection | **Diagnostic / experimental** | Summarizes chromosome-level topology but does not currently restrict optimizer states. |
| Bridge/articulation/2-core/biconnected decomposition of structural graphs | **Diagnostic / experimental** | Useful structural descriptors. They are not assumed to be objective-factorizing kernels. |
| Fundamental cycle basis | **Diagnostic / experimental** | Describes independent graph cycles without enumerating all simple cycles; not used in objective or branching logic. |
| Derived interchromosomal ordering constraints | **Infrastructure, not current hot path** | Implemented conservatively; current benchmark inputs do not yet provide enough genuine precedence/contiguity information to make this a major restriction. |
| Benchmark progress/checkpoint/resume | **Active validation infrastructure** | Long runs report per-case START/DONE state, checkpoint every completed case, and can resume without repeating finished cases. |
| Spectral ordering | **Deferred** | Not used by the current solver. |
| Braid-inspired moves / braid word metrics | **Deferred** | Not used by the current solver. |
| Tree-constrained ancestral event inference | **Deferred separate problem** | Downstream of stable extant-layout optimization. |

## Redundancy criterion

A stage is not called redundant because one easy benchmark did not require it.
The appropriate test is controlled ablation:

- hold the biological simulation and presentation seeds fixed;
- remove or bypass the candidate stage;
- compare correctness and proof status;
- compare lower-bound strength, residual state counts, nodes/transitions, and
  runtime;
- determine whether the stage provides unique audit or explanatory value.

If a stage supplies neither solver benefit nor unique interpretation, it should
remain diagnostic or be removed from the claimed core.

## Current graph-theory status

Graph theory remains genuinely active through the hypergraph incidence
representation and exact connected-component factorization.

The project now also has a second graph object: the **residual decision/factor
graph**. This graph is not a replacement for the biological hypergraph. It is a
derived representation of the remaining computational uncertainty after legal
reductions.

That distinction is important:

    biological hypergraph
        -> preserves observed n-ary homology evidence

    residual factor graph
        -> describes which unresolved legal decisions still interact in the
           crossing objective

The next graph-theoretic solver question is therefore not whether the raw
biological graph has an articulation point. It is whether the residual
objective graph has exact disconnected factors or sufficiently small
separators to permit proof-preserving decomposition.

## Parallelism status

Stage 18 introduces exact process-level parallelism across independent
incidence components in monotone branch-and-bound. The result is recombined
deterministically by adding component bounds and restoring component layouts in
canonical order.

The intended later hierarchy is:

    independent benchmark cases
      -> independent incidence components
      -> disconnected residual factor components
      -> separator-conditioned residual pieces
      -> parallel branch-and-bound frontier

Each layer must preserve the same global optimality proof: every unresolved
subproblem must be exhausted or bounded above the final incumbent before a
global optimum is declared.
