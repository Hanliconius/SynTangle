# Monotone solver pipeline

`optimize_auto` now uses one component dispatcher rather than the residual →
layer DP → branch-and-bound whole-fixture exception chain. Explicit standalone
exact APIs remain available, but `optimize_auto` does not restart through them.

## What survives a transition

| Reduction | Representation retained | Scope |
|---|---|---|
| Disconnected incidence components | One shared component map, solved once each; immutable results assembled at the end | Whole fixture |
| Hard orientation equations | Shared GF(2) basis passed to heuristic, factor builder and implicit search | Component |
| Factor-table invariance | Reduced table scopes consumed directly by continuation | That factor; never interpreted as a globally forced variable |
| Exact variable elimination | Reduced message tables plus live reconstruction records on the recursion stack | Current residual piece and branch |
| Separator conditioning | Conditioned tables and explicit condition context | Current branch only |
| Solved disconnected residual pieces | Their costs, bounds and assignments remain in the caller | Current branch |
| Bound-pruned branches | Removed from the continuation frontier; improving upper bounds only tighten pruning | Current piece/branch |

Factor materialization is preflighted for **all** tables before any scope
reduction. Oversized permutation domains or factor materialization therefore
route to implicit branch-and-bound before table eliminations begin, using the
same component map and orientation basis. Once factor reductions begin, cap
exhaustion is handled *inside* the current recursion by bounded search of those
exact reduced factors. There is no exception handoff back to the original
variables. Unexpected later cap exceptions fail loudly rather than silently
restarting.

Every factor-continuation handoff records active variables/scopes, retained
elimination records, retained conditions, factor-scope reductions, bounds and
whether it is unresolved. Runtime assertions forbid eliminated or conditioned
variables from re-entering that handoff. Conditional records are popped after
that branch, so they cannot exclude choices in another branch.

The continuation lower bound sums exact minima of each table consistent with
its partial assignment. It may be loose, but is admissible. Bounds add across
independent pieces; separator branches use the minimum of branch lower bounds.
Back-substitution reconstructs a complete legal layout, independently scored
against the unchanged biological evidence. Equality pruning retains one optimum;
this pipeline does not enumerate all tied optimal layouts. Node-budget exhaustion
leaves unresolved space with a bound, not a claimed exclusion or proof.

## First feasible layout

A dependency-free median homology-position sweep provides a GENESPACE-style
whole-chromosome order seed. This is our implementation of that general idea,
not a call to native GENESPACE and not an assertion of equivalent behavior.
Legal group flips and adjacent swaps refine it. Exact conditional order subset
DP in this preliminary pass is limited to 12 chromosomes per species/component;
large components retain all orders in the exact/bounded search. The implicit search also limits conditional order DP to 12 chromosomes,
using an admissible independent-pair relaxation above that size. This avoids
re-entering the same exponential DP in its lower-bound calculations.
The preliminary
pass is capped at 20 improving iterations per restart. A heuristic decision
is never treated as a proof that another order is unnecessary.

This changes search strategy; identical layouts or evaluation counts are not
expected relative to the previous cached scorer release. Compare objective,
valid bounds, time, and proof status instead.

## Short Pegasus comparison

Use an isolated worktree and run `bash validation/benchmark/submit_pipeline_benchmark.sh`.
The nine unchanged stress inputs run as a Slurm array. Each task compares the
previous cached release with this pipeline, at 120 seconds per method, one
restart, five improving heuristic steps, and 1,000 bounded-search nodes per
component. The complete auto solver is timed, not only local scoring.

`paired.json` retains completed results, bounds, reduction audits and independently
validated incumbent checkpoints for timeouts. Prior native GENESPACE and
GENESPACE-plus-our-flips scores are read from the original stress comparison;
they are not recomputed. Original inputs and ongoing result files are read-only.
A timeout is never counted as an optimization win or a proven optimum. These
short runs establish whether useful incumbents/results arrive sooner; they do
not replace the longer quality/proof experiments.
