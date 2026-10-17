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
| Recursive residual factor elimination | **Active primary exact solver machinery** | Builds exact crossing-factor tables, removes provably invariant scope variables, eliminates objective-neutral/leaf variables, splits disconnected residual pieces, conditions on small articulation variables, and uses min-fill elimination only on the irreducible remainder. Independent incidence components can run in separate worker processes. |
| Exact species-layer dynamic programming | **Active exact fallback** | Historical exact solver retained as an independent formulation and fallback when residual intermediate tables exceed configured limits. |
| One-layer subset DP | **Active solver machinery** | Used by local search and branch-and-bound for conditional chromosome-order optimization and bounds. |
| Constraint-aware local search | **Active supporting machinery** | Produces legal feasible incumbents and therefore useful upper bounds. |
| Monotone component branch-and-bound | **Active solver machinery** | Handles cases beyond exact layer-DP limits; reports explicit LB/UB/gap and can parallelize independent incidence components. |
| Residual variable/factor graph | **Active solver representation** | Represents the unresolved legal order/orientation decisions and the adjacent-layer crossing factors that actually couple them. The graph is rebuilt/reduced through exact table operations rather than used only as a diagnostic. |
| Residual disconnected-piece solving | **Active solver machinery** | Disconnected residual factor pieces are solved recursively and their optima added exactly. Newly disconnected pieces created by elimination/conditioning are split again rather than rejoined. |
| Objective-neutral / leaf elimination | **Active solver machinery** | Conservative factor scopes are contracted when exact factor tables prove a variable irrelevant; variables with residual primal degree at most one are peeled immediately by exact min-sum reduction, including endpoints carrying unary message factors. |
| Small articulation separators | **Active solver machinery** | Small-domain articulation variables are conditioned exactly, exposing independent subproblems that are solved and recombined for each legal separator value. |
| Min-fill variable elimination / treewidth-aware reduction | **Active exact fallback inside residual solver** | When no cheaper split remains, an exact min-fill elimination step contracts one variable. Intermediate table/work caps prevent uncontrolled blow-up. |
| Larger separator / tree decomposition | **Experimental next step** | Multi-variable separators and richer tree decompositions remain candidates when single articulation conditioning and min-fill elimination are insufficient. |
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

The active graph-theoretic solver question is therefore not whether the raw
biological graph has an articulation point. It is whether the residual
objective graph has exact disconnected factors or sufficiently small
separators to permit proof-preserving decomposition. Stage 20 now acts on
those properties directly: exact factor tables can shrink conservative scopes,
leaf variables are eliminated, disconnected pieces are solved separately, and
small articulation variables are conditioned before generic elimination.

## Parallelism status

Independent incidence components can now be dispatched to worker processes by
both the primary residual exact solver and monotone branch-and-bound. Inside
each incidence component, the residual exact solver recursively exposes finer
objective decomposition.

The active hierarchy is:

    independent benchmark cases
      -> independent incidence components
      -> exact factor-scope contraction
      -> disconnected residual factor components
      -> objective-neutral / leaf elimination
      -> articulation-conditioned residual pieces
      -> min-fill exact elimination
      -> legacy exact DP or monotone B&B only if configured residual caps fail

Each exact layer preserves the same global optimum. Bounded branch-and-bound
still requires every unresolved subproblem to be exhausted or bounded above
the final incumbent before a global optimum is declared.
