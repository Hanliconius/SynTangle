# Pegasus benchmark history and solver decisions

Updated 2026-10-06. This is the development record for the
`pegasus-genespace-pilot` branch. The joint and hybrid solvers are experimental;
these results have not switched the default production pipeline.

## Evidence and interpretation

C is the canonical count of adjacent-species homology-link crossings under legal
whole-chromosome ordering/reversal. Internal order and observed homology remain
fixed. A layout is an upper bound; a matching valid global lower bound establishes
optimality under this model. A drawing optimum does not establish an ancestral
history or biological accuracy. Zero crossings is optimal by nonnegativity.

Native GENESPACE and **GENESPACE plus flips** are distinct comparisons. The latter
is our external assistance routine, not native GENESPACE functionality. The
native layout uses a fixed reference row; SynTangle may optimize all rows. The
comparison concerns the shared crossing objective, not every GENESPACE feature.

The record below comes from user-supplied Pegasus outputs. Compact collector
results are archived where available. Full generated bundles, layout JSON,
solver logs and fingerprints remain on Pegasus until explicitly exported. Do
not describe the archived stdout as an independently checked proof certificate.

## Experimental progression

| Experiment | Known Pegasus jobs | Observation | Decision |
|---|---|---|---|
| Paired smoke | 138371221 | Three presentations, initial C=390/912/1064, all C*=0 in 0.22–0.33 s; paired biological evidence unchanged | Input/presentation invariance smoke passed |
| GENESPACE environment | 138371233 | GENESPACE 1.3.1 available | Proceed with native comparison |
| Layout pilot and largest-case retry | 138371248; collector 138371249; retry collector 138372078 | Ten complete cases: five ties and five SynTangle improvements versus external flip assistance | Native GENESPACE does not always minimize this objective; assistance remains a separate baseline |
| Coupled stress baseline | 138372260; collector 138372261 | 6-species triplet reached C*=341 in about 720–746 s; larger tasks ran >2 h and were cancelled | Profile bottlenecks rather than increase resources |
| Scoring cache | 138372623; longer probe 138372720 | 4.8–5.0x local-search speedup at 6 species, 7.0–7.4x at 8; identical search in completed comparisons; 10-species cases still timed out | Keep cache; investigate a second bottleneck |
| Large-case diagnosis | 138373188 | Stack samples inside conditional subset-order DP | Cap exponential conditional DP; preserve unresolved choices |
| Monotone pipeline | 138373390 | Faster small cases and better large-case incumbents, but larger runs still externally timed out | Keep shared preprocessing and reduced-factor continuation |
| Cached branch machinery | Job ID not retained in pasted evidence | Identical result/search; ~1.34x at 8 species and ~3.95–3.97x at 10 | Keep exact caches |
| Search budgets/deadlines | 138373671; 138373729 | Cooperative deadlines returned layouts/bounds instead of losing results; more search improved large layouts, but proof gaps stayed large | Keep checkpoints and explicit bounds; avoid equating completion with proof |
| Orientation clusters | 138383947; implementation f9a9b62 | No material bound/score gain; fewer branch nodes on largest cases | Leave default off |
| Conditional-order probe | 138384046; implementation f7228c8 | Substantial conflicts between the two neighbors' preferred shared order; sampled triangle lifts zero | Target shared ordering across layers; do not treat conditional bounds as global |
| Coupled factor buckets | Implementation a3a57dc; job ID not retained | 6-species lower bound 276→298; 8-species +1; largest branch work ~1571→122 without better layouts | Park bucket refinement; cost outweighed pruning |
| Global methods | Array 138387348, setup 138387347, collector 138387349; implementation 5618814 | Joint MILP established reported C*=341 at 6 species and C*=1258 at 8; neighborhoods best on largest cases | Prioritize joint formulation and hybrid neighborhoods |
| Hybrid optimization and strict-improvement audit | Implementation a8c0bde; current job IDs not supplied | Submitted by user; results and audit output pending | No outcome recorded yet |

Earlier pilot task 138371248_9 failed with exit 124 after 25:01 and ~435 MB RSS,
with an interrupted-system-call message. Accounting initially showed a stale
PENDING record; `sacct --duplicates` exposed the failure. The timed-out run is
not classified as a win. Its retry completed: 8 species/31 chromosomes/zero
planned events, native GS C=6996, assisted C=1320, SynTangle C*=0 in 0.4508 s.

## Stress inputs held fixed

Chromosome counts below are ancestral counts per genome; fusion/fission can
change extant counts. Events are per lineage step. These rungs jointly increase
species, chromosomes and event burden, so they do not isolate the separate
causes of scaling. Each rung has three presentations of the same biology.

| Rung | Species | Ancestral chromosomes | Events/step | Anchors/chromosome | Biology seed |
|---|---:|---:|---:|---:|---:|
| Small | 6 | 20 | 5 | 16 | 8101 |
| Medium | 8 | 30 | 8 | 18 | 8201 |
| Large | 10 | 40 | 12 | 20 | 8301 |

The Pegasus input root is `SynTangle_pegasus_test/local_results/genespace_stress`.
Use its existing manifest and case directories; do not regenerate biology for
method comparisons.

## Global-method result: the formulation change

The [archived collector summary](../validation/benchmark/results/2026-10-06-global-methods/summary.md)
contains all 108 rows, with [machine-readable records](../validation/benchmark/results/2026-10-06-global-methods/records.json)
and [provenance](../validation/benchmark/results/2026-10-06-global-methods/metadata.json).
There are 106 complete records and two SDP wall cutoffs. Complete means a result
was returned, not necessarily that an optimum was proven.

At a 150-second allowance:

| Case/presentation | Starting C | Pipeline C | Joint MILP C | MILP lower | MILP wall s | Neighborhood C | SDP-rounded C |
|---|---:|---:|---:|---:|---:|---:|---:|
| 6sp mild | 341 | 341 | 341 | 341 | 1.02 | 341 | 341 |
| 6sp strong | 341 | 341 | 341 | 341 | 1.17 | 341 | 341 |
| 6sp random | 341 | 341 | 341 | 341 | 1.02 | 341 | 341 |
| 8sp mild | 9440 | 9440 | 1258 | 1258 | 2.72 | 5306 | 2333 |
| 8sp strong | 5038 | 5038 | 1258 | 1258 | 2.92 | 4991 | 2286 |
| 8sp random | 7931 | 7931 | 1258 | 1258 | 3.19 | 4490 | 2623 |
| 10sp mild | 49846 | 49619 | 38084 | 4343 | 27.61 | 25700 | 48424 |
| 10sp strong | 59043 | 57913 | 40704 | 4343 | 27.40 | 23272 | 57043 |
| 10sp random | 49336 | 48672 | 34574 | 4343 | 27.96 | 21804 | 47336 |

The same reported medium optimum across all three presentations is an important
invariance result. Relative to the pipeline layouts, the medium MILP removes
75–87% of crossings; large-case neighborhoods remove 48–60%. Large-case optima
remain unresolved. These observations support a formulation change, not a
claim that the exact solver scales to all biological datasets.

The original MILP divided time equally between remaining components, often
visiting a hard component before trivial ones. Its 600-second large runs finished
in ~93–96 seconds with the same scores as the 150-second runs. They do not test
a full ten minutes of global search. Neighborhoods did use ~600 seconds but
returned the same scores as at 150 seconds. The hybrid experiment corrects time
allocation and stops repeating unchanged neighborhood sweeps.

SDP used a full Shor relaxation only on small decision models and a shared-moment
block relaxation on larger ones. Its numerical objectives are not certified
bounds. Weak block-SDP performance does not rule out the stronger published
multi-level SDP formulation.

## Current hybrid experiment

See [HYBRID_METHODS.md](../validation/benchmark/HYBRID_METHODS.md). Four variants
use the same direct HiGHS version: corrected scheduling without a MIP start,
with a MIP start, with three-species neighborhoods, and with adaptive three-to-five
species neighborhoods. Each is tested on all nine existing cases at 30/150/600 s,
one CPU and 16G per task, without an array throttle. Starting layouts are frozen
from the latest best validated saved outputs, including the global-method run.

Zero-crossing independent components are retained as proven optimal. Other
components reuse graph/GF(2) reductions, crossing factors and sparse model
matrices. Conditional neighborhood restrictions do not become global exclusions.
Each component enters global exact search once; no branch-tree restart occurs
inside the hybrid. This does not implement import of an old factor-elimination
or branch-search frontier.

Setup also re-scores previous MILP results and tests whether C <= reported_C - 1
is feasible for each positive-score component of a reported optimum. Counterexamples
abort setup; time-limited audit attempts remain unresolved. The audit is another
numerical solver check, not an independently checked rational proof file.

Audit/collector results will be catalogued after they are supplied. No job ID,
new optimum, or claimed speedup is inferred from a submission instruction alone.

## Promotion and remaining validation

Promote a strategy only after legal-state/objective checks, reported-bound audit,
matched-budget comparison, presentation invariance and visual reconstruction
checks. The hybrid implementation currently has 111 passing unit tests, including
exhaustive objective equivalence, transitivity, hinted solver reconstruction,
strict-improvement checks and removal of proven-zero components.

Still outstanding: more independent biology seeds, the orthogonal scaling
profile, real synteny inputs, visual audits of the newly improved layouts, and
an explicit residual-domain/reconstruction contract before any backend handoff.
