# Solver efficiency ablation

This measures contributions to the joint solver's improvement; it does not
promote a new production default. All big runs are submitted to Pegasus.

## Controlled comparisons

| Variant | Global backend | Allocation among unresolved components | Feasible MIP start | Neighborhoods |
|---|---|---|---|---|
| equal_nohint | direct HiGHS | equal remaining share | no | none |
| weighted_nohint | direct HiGHS | decision/product-count weighted | no | none |
| weighted_hint | direct HiGHS | weighted | yes | none |
| weighted_three | direct HiGHS | weighted | yes | three adjacent species |
| weighted_adaptive | direct HiGHS | weighted | yes | three, then five species |
| weighted_scipy | SciPy MILP | weighted | no | none |

Equal vs weighted isolates allocation **within the new shared model**. Both
remove proven-zero components and visit smaller unresolved models first. It
does not recreate the old scheduler's ordering or time wasted on zero components.
SciPy vs direct HiGHS also includes API/model conversion differences, not just
a solver-version difference. Caches are retained in all variants; their benefit
was already measured with identical search in the earlier paired experiments.

Two starting states are frozen per fixture:

- `public`: public chromosome order with the hard GF(2) base orientation
  assignment applied, followed by component canonicalization. No saved optimized
  layout is used. This is a solver cold start, not a complete GENESPACE or
  default-pipeline end-to-end comparison.
- `pre_hybrid`: the historical frozen input from the earliest available hybrid
  run, re-scored and validated. It deliberately does not select today's best
  known optimum. Override the run directory with `SCT_WARM_RUN` if necessary.

All nine existing presentation fixtures receive all six variants and both starts
at 150 seconds, twice: 216 tasks. The three largest receive the same comparisons
at 600 seconds once: 36 tasks. Total 252, one CPU and 16G each, no array throttle.
Repeats use identical inputs and options to measure runtime variability; they
are not independent biological replicates. Single 600-second timings cannot
establish a reliable timing ranking.

## Measurements and proof scope

Every final result includes graph-model preparation, sparse matrix construction,
neighborhood time/subsolves, global call time, backend solver time, global branch
nodes, peak worker RSS, versions and incumbent improvement checkpoints.
Per-component details include problem size, final bounds and neighborhood work.
Stage timings are inclusive and overlap; do not sum every timing column.
Worker subprocess wall time includes interpreter/import/input overhead. Frozen
start preparation and input fingerprinting happen in setup and are separate.
SciPy supplies returned incumbents only; direct HiGHS also checkpoints improving
callbacks. The trajectory therefore cannot fairly compare time to every
intermediate incumbent across the two APIs. Bound progress is recorded at
component completion, not sampled during the global search.

All runs retain the same incidence/GF(2) reductions, sparse crossing factors and
hard constraints. Temporary neighborhood restrictions never become global
exclusions. Each component starts global search once; changing variants means
independent experimental runs, not restarting a pruned search frontier.

A separate nine-task audit array re-scores historical global and hybrid results,
then asks whether a strictly better layout is feasible for each positive-score
component of a reported optimum, with 120 seconds per component. Zero-crossing
components are checked by nonnegativity. Counterexamples fail the audit;
unresolved attempts and missing audit records remain explicit. This is another
numerical solver check, not an independently checked rational proof certificate.
Audit runs concurrently with efficiency runs; collector waits for both arrays.

## Submit and inspect

    bash validation/benchmark/submit_efficiency_benchmark.sh

The launcher snapshots source/runner code, freezes starts and input file hashes,
and records commit and job IDs in a unique `local_results/efficiency_*` directory.
The setup installs the same pinned dependencies as the hybrid experiment.

    cat logs/st_eff_setup.*.{out,err}
    cat logs/st_eff_audit.*.{out,err}
    cat logs/st_eff_collect.*.{out,err}

The collector emits an overview plus per-case results. Machine-readable records
retain all stage/component details and audit checks. Completion alone is not
optimality; only matching valid global lower/upper bounds support that claim.
Missing, failed and cutoff tasks are never wins.

## Interpretation

Compare equal_nohint with weighted_nohint for allocation, weighted_nohint with
weighted_hint for hints, weighted_hint with the two neighborhood variants for
heuristic work, and weighted_nohint with weighted_scipy for the API/backend.
Compare the two frozen starts within each matched case/variant/budget. Record
layout quality at the cutoff separately from time to matching bounds. Repeat
results and stage costs will guide defaults; avoid selecting by one wall time.
