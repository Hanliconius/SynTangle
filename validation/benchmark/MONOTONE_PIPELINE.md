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

## Diagnose unresolved large cases

`bash validation/benchmark/submit_pipeline_diagnostic.sh` submits a read-only
saved-audit job followed by three 90-second sampled solves: 8sp mild, 10sp mild,
and 10sp strong. It uses the same five-step/one-restart/1,000-node settings as
the short comparison. Solver code and decisions are unchanged; wrappers record
call counts, inclusive/exclusive timing, orientation pruning/forcing, and
samples of elimination, scope reduction and conditioned continuation context.
Faulthandler samples the live Python stack every 30 seconds. A controlled
90-second interruption saves the counters and active stage even without a
final solver result. These instrumented times are diagnostic, not speedup
measurements. Missing final audits from earlier timeout jobs are explicitly
reported as unknown, never taken as evidence of successful pruning.

## Cooperative deadlines and live incumbents

`optimize_auto(..., time_limit_seconds=150, component_workers=1)` now
returns an ordinary bounded result when its cooperative deadline is reached.
Branch-and-bound checks the deadline in orientation reduction, at accepted
nodes, and during permutation scans (including candidates rejected by bounds).
An incomplete ordering search retains its relaxed lower bound; the orientation
frontier retains its open bounds. Completed components remain completed.
Unstarted components keep the legal heuristic incumbent and lower bound zero.
Deadline exhaustion is never recorded as a proof-based exclusion.

Residual elimination can switch, at an elimination boundary, to bounded
continuation on its current reduced factors. An expired deadline prevents
further continuation expansion; eliminated variables and conditioning scopes
are reconstructed through the existing stack. This does not restart the search.
Preparation, factor materialization and individual scoring/DP operations are
not asynchronously interrupted: the deadline is cooperative, not a hard wall
clock guarantee. An external emergency cutoff remains useful as a backstop.
Deadline-controlled multiprocessing is explicitly unsupported at present.

Each strictly improved complete feasible component state discovered inside
implicit ordering search is merged with the best retained states of all other
components, canonically rescored, and sent to the existing progress callback.
This saves improvements even if the external emergency cutoff eventually
terminates the worker. Residual improvements are published when that component
returns; individual partial factor assignments are not exported as layouts.

`submit_search_budget_benchmark.sh` now launches all 18 tasks concurrently:
500 and 5,000 nodes per component on each unchanged stress input, with a
150-second cooperative deadline and 210-second emergency process cutoff.
Completed records include layout, valid bounds, optimality gap and a deadline
flag. Emergency cutoffs retain the last exported feasible incumbent and make
no claim about a final bound.

## Saved incumbents and coupled orientation bounds

`starting_state` accepts a complete saved `LayoutState`. Before use it must
contain every species and chromosome exactly once, have only signs +1/-1,
and satisfy the shared hard GF(2) orientation basis. Its canonical score is
compared with the median seed; the better seed starts local search. Final
upper bounds cannot exceed the retained incumbent. A saved feasible layout
is an upper bound, not a proof that any other ordering is infeasible.

`bound_cluster_size=6` enables an experimental stronger orientation lower
bound (default 0 preserves the independent-cell bound). Groups are partitioned
deterministically into clusters of at most six. Each internal anchor-cell
factor is assigned once: if its entire scope lies inside a cluster, its cost
is summed with the other internal factors and minimized jointly over compatible
cluster bits. Factors spanning clusters retain their independent relaxation.
Complete and partial cluster tables are precomputed, preserving shared-group
correlations and strict floating-point coordinate ties. Summing these independent
cluster/cross-cluster minima remains admissible and cannot be weaker than the
old per-cell relaxation for the same partial assignment. It does not assert
chromosome-order consistency or discard branches without a bound proof.
The 6-group table has at most 729 partial assignments. The existing memory
request is unchanged. Cluster-controlled solves currently require one worker.

`submit_bound_pair_benchmark.sh` launches nine cases concurrently. Each case
selects and rescoring-validates the best saved layout from prior benchmark
folders, then runs cluster sizes 0 and 6 from that identical layout, at 5,000
nodes per component and a 150-second cooperative deadline. This compares
algorithmic search changes, so identical search counters are not expected.
Invalid saved layouts or inconsistent recorded scores fail explicitly.
