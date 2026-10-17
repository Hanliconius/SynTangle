# Shared-model hybrid optimization experiment

The global-method comparison showed a step change for joint MILP: medium stress
cases reached the same reported optimum C=1258 across three presentations in
about three seconds, while graph neighborhoods produced the best large-case
layouts. This experiment tests their combination on unchanged input evidence.
The default production pipeline is unchanged.

## Four controlled variants

All global solves now use the same pinned direct HiGHS backend (highspy 1.15.1),
one CPU, and a common frozen best saved incumbent for each fixture:

1. `milp_reclaim`: correct component scheduling; objective cutoff and retained
   incumbent, without passing a feasible starting solution to the solver.
2. `milp_hint`: same, with a complete feasible MIP start including product variables.
3. `hybrid_3`: same feasible start, preceded by adjacent three-species neighborhoods
   on components with more than 128 decisions.
4. `hybrid_adaptive`: expands from three to five adjacent species after an
   unchanged sweep. Neighborhoods jointly release chromosome precedences and
   all hard flip groups touching those species.

The direct backend differs from the HiGHS version bundled in SciPy, so comparisons
with the previous experiment include that backend change. Within this new array,
the global backend is identical across all variants. Neighborhood subsolves use
the cached SciPy MILP representation, matching the earlier neighborhood method.

## Shared work and exclusion scope

The graph components, GF(2) basis, aggregated crossing factors and immutable
sparse MILP matrix are constructed once per component and reused. Neighborhood
bounds/fixings have conditional scope and are never turned into global exclusions.
A sweep that cannot improve stops repeating; this is a heuristic stopping rule,
not a claim that the global region has been excluded. The adaptive variant then
expands its neighborhood.

Every component gets exactly one global exact search after its heuristic stage.
There is no exact-search restart, no transfer to another backend after branching,
and no loss of exclusions from that global solver's presolve/branch tree during
its run. Its best incumbent and lower bound are exported when it returns. This
is still a fresh experiment from shared graph preprocessing, not an importer for
an old factor-elimination or branch-search frontier.

## Time use

Independent components already at zero crossings are retained as proven optimal
and removed before model construction; their layouts remain reconstructable.
Remaining components are scheduled from smaller to larger models. Remaining time is weighted
by decision/product count and recalculated after each completed component. The
last unresolved component receives all remaining time. This fixes the old policy
that gave the large component only its equal share before visiting many trivial
components, then finished with unused time.

For hybrid variants, neighborhoods use at most 25% of their component allocation
and at most 60 seconds, with each conditional solve capped at five seconds. The
remaining component allocation goes to one global solve. Model preparation and
matrix construction count toward elapsed work. Small models skip neighborhoods.

## Proof audit and checkpoints

The setup job re-scores all previous complete MILP results against current input.
For each case with a reported optimum it checks a representative optimal layout,
sums component scores, and solves each positive-score component with the added
constraint C <= reported_C - 1. Infeasibility validates the reported optimum;
a feasible counterexample aborts setup; an audit timeout is explicitly unresolved.
Zero-score components use the mathematical nonnegative objective bound.

This audit uses another HiGHS invocation, not an independently checked rational
proof file. It tests the strict-improvement feasibility formulation and canonical
reconstruction. It runs only as verification, never as a solver fallback.

Direct HiGHS exports improving legal layouts during global search via callbacks;
the worker canonically scores each one and writes atomic checkpoints. Final
records include raw dual bound, objective constant, component bounds, solver
status/version, preparation and stage times, and whole-chromosome changes.

## Pegasus benchmark

Nine existing cases x four variants x 30/150/600 second limits = 108 independent
tasks, submitted without a throttle. Each uses one CPU and 16G RAM. Dependencies
are installed in an isolated venv by a compute-node setup job. This new array
includes the best saved layouts from `global_methods_*`, so it asks for further
improvement, rather than starting over from the older inferior incumbents.
All methods/budgets freeze exactly the same starting layout per case.

    bash validation/benchmark/submit_hybrid_methods_benchmark.sh

Audit/start preparation log:

    cat logs/st_hybrid_setup.*.{out,err}

Final comparison:

    cat logs/st_hybrid_collect.*.{out,err}

Missing/failed/cutoff results are not wins. Global and conditional bounds remain
separate. SDP is parked for this comparison because its larger-model cost did
not match the benefit observed from MILP and neighborhood search.
