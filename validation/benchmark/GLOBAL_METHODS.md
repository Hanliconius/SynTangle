# Joint formulation experiment

Four independent experiments on the existing nine stress fixtures, each at
30, 150 and 600 seconds. All 108 tasks use one CPU and 16G RAM and the same
frozen, canonically validated best saved layout for each case. No simulation
is regenerated. GENESPACE/native and external flip-helper scores are retained
in each task record. The pipeline uses a generous node ceiling so time, rather
than the previous 5,000-node cap, normally limits it.

## Methods

- `pipeline`: current default solver, experimental cluster bounds disabled.
- `milp`: joint binary chromosome precedences and propagated GF(2) flip groups.
  Transitivity constraints define total chromosome orders. Aggregated crossing
  tables become an exactly equivalent quadratic objective; continuous product
  variables with McCormick constraints linearize it exactly for binary decisions.
  SciPy/HiGHS solves the MILP. The saved layout supplies an objective cutoff and
  retained feasible solution, not a solver solution hint. Dual bounds are reported
  with a conservative integer rounding tolerance; these are ordinary numerical
  solver certificates, not independently checked rational proof files.
- `lns`: the same exact model restricted temporarily to three adjacent species
  and all hard flip groups touching them. Each subsolve gets at most 10 seconds;
  neighborhoods rotate through the species. Outside decisions are fixed to the
  current legal layout. Restricted bounds are never exported as global bounds.
- `sdp`: Shor semidefinite relaxation of the same binary objective with shared
  first moments, order-transitivity constraints and McCormick inequalities. For
  components with <=128 decisions it uses one full PSD matrix. Larger components
  use connected blocks of <=32 decisions and 3x3 PSD constraints for cross-block
  objective pairs. This is an implementable *block relaxation*, not a reproduction
  of Chimani et al.'s stronger published formulation. All objective factors remain
  present. SCS numerical objectives/statuses are diagnostic estimates only, never
  certified bounds or pruning decisions. Deterministic/fixed-seed randomized
  projections onto legal whole-chromosome orders are canonically scored, retaining
  the saved layout when rounding does not improve it.

## Reductions and evidence

Each method branches at the shared preprocessing boundary, before any existing
exact-search exclusions. Disconnected incidence components and hard orientation
relations are reused via the same implementation. Zero-cost pair decisions are
included when needed to enforce total orders. There are no permutation domains.
Only whole chromosomes move or reverse. Every reconstructed state is validated
against hard constraints, and MILP objective values are checked against canonical
crossing scores. Output states and component group mappings are recorded.

This experiment does **not** promise to import an existing reduced-factor search
frontier or old branch exclusions. It compares independent solves from the same
preprocessed problem; it is not a fallback that restarts an earlier search. A
future production integration must explicitly translate residual domains,
eliminations, objective constants and reconstruction ledgers before any handoff.
Neighborhood fixings are conditional and never permanent global reductions.

## Timing and failure handling

Each method starts afresh for each budget: longer runs do not import outcomes
from shorter runs. Report end-to-end worker wall time, model preparation time,
score, global bounds where available, SDP compilation time and failure status.
Preparation/compilation counts toward the allotted work; the parent allows 20
seconds of finish/reporting overhead before an emergency process cutoff. Saved
legal checkpoints survive that cutoff. Internal MILP incumbents appear only when
a subsolve returns; their intermediate improvements cannot be recovered after
an emergency kill. Missing or failed results are not wins. The collector counts
records against the expected 108 tasks and includes statuses, scores and bounds.

A dependency/preparation Slurm job creates a separate venv within the unique
results directory, installs pinned dependencies, and freezes all nine inputs'
starting layouts. The array depends on its success. No existing micromamba env
is modified; installations and input validation run on a compute node. All tasks
are submitted together; cluster availability determines actual concurrency.

Run:

    bash validation/benchmark/submit_global_methods_benchmark.sh

Then inspect:

    cat logs/st_methods_setup.*.{out,err}
    cat logs/st_methods_collect.*.{out,err}

Research basis: quadratic linear ordering (Buchheim, Wiegele & Zheng),
PACE 2024 MPPEG one-sided crossing branch-and-cut, local branching
(Fischetti & Lodi), and multi-level SDP (Chimani et al.). Performance of those
methods on their problems does not establish performance on these fixtures.
