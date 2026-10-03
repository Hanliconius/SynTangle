# Monotone search-space reduction

SynTangle's optimization pipeline follows normative rule R17: once a
decision state has been proved impossible, equivalent, forced, dominated by a
valid bound, or unable to improve the incumbent objective, later stages must
not recreate it.

The distinction is between **biological evidence** and **search variables**.

- Full anchor/block evidence is retained for exact scoring, reconstruction,
  and audit.
- Only the unresolved legal display decisions are reduced.

The intended flow is:

```text
full extant evidence
        |
        v
independent chromosome-homology components
        |
        v
hard orientation propagation (GF(2))
        |
        v
compact free-flip basis
        |
        v
safe lower bounds + incumbent
        |
        v
force/prune orientation alternatives
        |
        v
branch only on residual orientation groups
        |
        v
fixed-orientation order search
        |
        v
force/prune order alternatives
        |
        v
proven optimum OR explicit lower/upper bound gap
```

## Stage 15 implementation

The branch-and-bound solver no longer converts the compact orientation basis
into every complete orientation assignment before order optimization.

Instead it keeps one ternary residual vector per free orientation group:

- `None`: unresolved;
- `0`: retain the basis orientation for the group;
- `1`: flip the complete group.

A safe lower bound is computed by relaxing global chromosome-order
consistency. Link-pair crossings are partitioned into disjoint classes:

1. links on the same chromosome pair contribute their minimum possible
   within-chromosome crossing count under unresolved legal orientations;
2. links sharing the left chromosome but ending on different right
   chromosomes contribute the better of the two pairwise right-chromosome
   precedences;
3. the symmetric shared-right-chromosome case is treated analogously;
4. links whose two chromosome endpoints are both different contribute zero in
   the relaxation.

Because pairwise precedence choices are allowed to disagree with one another
in the relaxation, this bound can be optimistic but cannot exceed the true
constrained optimum.

At every residual orientation node, SynTangle compares both values of each
unresolved flip group against the current incumbent. If one alternative has a
lower bound at least as large as the incumbent, that alternative is removed
and the other value becomes forced. Propagation repeats to a fixed point.

When an incumbent improves, queued nodes are propagated again against the
new, tighter upper bound before further branching.

## Reuse rather than re-solving

Stage 15 caches:

- lower bounds for canonical partial orientation assignments;
- exact one-sided order-DP subproblems keyed by orientation, target layer,
  neighbor layer, and fixed neighbor order;
- canonical residual orientation states already placed in the search queue.

This avoids repeating calculations when different search paths expose the
same residual subproblem.

## Decomposition

Disconnected chromosome-homology components remain independently solvable
under R7. Further dynamic splitting is allowed only when independence of the
**residual objective factorization** is proved. Stage 15 does not infer
independence merely because orientation variables are absent from the same
local bound term: chromosome-order variables can still couple them.

This is intentionally conservative. A future factor-graph decomposition may
split residual kernels further, but it must prove that no remaining objective
term spans the proposed split.

## Diagnostics

Per component, branch-and-bound now reports:

- number of free orientation groups;
- number of complete orientation states represented implicitly;
- root relaxed lower bound;
- unresolved orientation groups after root propagation;
- orientation nodes actually evaluated;
- complete orientation leaves that reached order search;
- orientation branches pruned;
- groups fixed by bound propagation;
- order-search nodes evaluated;
- memo/cache hits;
- final lower bound, upper bound, gap, and proof status.

These diagnostics make the reduction process auditable: the solver reports
where combinatorial complexity was removed rather than only the final score.
