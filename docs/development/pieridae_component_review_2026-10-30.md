# Pieridae: component-first diagnostic and repair experiment

## What was inspected

SynTangle main `f094d75f881913373126b3c56f9bef0ca4329a72`, the
conditional/interface/backbone experiments in PRs 26–28, MerianCode main
`b85a0272e16404e166eec44107d4fb60e2b4eca9`, ColiasGenomics_2026 main and
`pieridae-supplied-topology-20261029`, and the supplied two- and 31-species PDFs.
The original Pegasus tables, solver logs and checkpoints are not available in
this review environment. Therefore the actual Z component and real-data runtime
improvement remain unverified. PDFs cannot establish graph connectivity or an
optimality certificate.

The 31-species PDF shows 14,170,379 crossings on the left and 14,889,733 on the
right. The right panel is explicitly an ordering-only backbone heuristic with
all links restored. The code optimizes a support/run-filtered graph and accepts
that proposal even if the full score increases. This is a different objective,
not failure of a completed full-evidence exact solve. The right result is 719,354
crossings worse on the stated full objective and must not replace the better
incumbent in a canonical run.

The two-species PDF compares an all-evidence bounded solve (682,299) with a
masked-evidence bounded solve rescored on all links (690,711). It is not an
input-versus-optimized comparison, nor evidence that either solve is globally
optimal. Both plots preserve an important visual observation: the largest
changes cluster at interfaces between chromosome organizations.

MerianCode's current `syntangle_plot.py` calls the public `syntangle run` entry
point. It does not implement its own chromosome optimizer. The public run uses
`run_hybrid`; the earlier Colias initial pipeline used `optimize_auto`, while
later refinements use `run_hybrid`. Those paths have different preparation and
incumbent-generation behaviour. Do not extrapolate this fix to every path.

## Confirmed implementation issues and distinctions

1. **Preparation can dominate search.** `CoupledCrossingBound` constructs factors
   by enumerating pairs of links at every interface. With m links this visits
   m(m-1)/2 pairs before MILP search. The old lower-bound kernels also explicitly
   enumerate anchor pairs. Solver time limits do not cap those loops.
2. **Easy structure is handled too late.** The hybrid path already canonicalizes
   components and skips components whose incumbent is zero, but does not have a
   dedicated exact presolve for a long strand with intrinsic crossings. It also
   builds all component models before solving the first one. The alternative
   pipeline runs whole-fixture local search before exact component dispatch.
3. **Connectivity must use full evidence and hard constraints.** One retained
   Z/autosome link can join their components. This might reflect ambiguity,
   duplication, input mistakes or biology; layout code cannot decide to delete
   it. Filtered-graph independence is not full-graph independence.
4. **Component order is not pixel alignment.** A common component order removes
   interleaving and cross-component inversions. Different chromosome lengths,
   absent components and row normalization can still shift its physical x
   position. An optional component-lane renderer could reserve a common span
   per component with explicit blank space. Do not alter genomic coordinates
   to make a chromosome appear aligned.
5. **Species adjacency is not evolutionary time.** A line between two extant
   rows is a comparison, not an ancestral transition. Nesting a Pieris-like
   group between two ancestral-like groups can display two boundaries without
   implying two independent historical reorganizations.

## Implemented experiment

Opt in with `syntangle run --component-presolve ...`, or the Pegasus runner
below. No normative rules or existing defaults are changed, except replacement
of two lower-bound counting kernels by equivalent faster counting algorithms.
No anchors are removed, sampled or moved internally.

Before general model construction, the hybrid path now:

- Builds components from the full chromosome–homology incidence graph, merging
  components linked by hard orientation constraints as before.
- Canonicalizes their common order across all rows.
- Solves a component with at most one chromosome per row by an exact two-state
  orientation dynamic program. Nonadjacent hard equations are deferred to the
  general solver. Empty rows remain empty; no new adjacency is introduced.
- Retains any other feasible zero-crossing component, using the universal zero
  lower bound. Both types are removed from later model construction.
- Builds unresolved component models one at a time, after these easy results
  have been checkpointed. The experiment uses lossless chromosome-pair-cell
  aggregation and a cooperative deadline in its cell-factor loop. If preparation
  stops, it retains the incumbent and reports zero as the unresolved component's
  bound; it never exports an unfinished factor model as a proof.

The general solver remains the existing full order-and-flip MILP with adaptive
neighborhoods. There is no silent heuristic replacement. This is a fresh run
from a feasible saved layout, not continuation of an earlier proof frontier.
The package's numerical bound closure is still distinct from an independently
checked exact certificate.

### Why component separation and the strand recurrence are valid

For disconnected decision components, restrict any full layout to each
component. Its internal crossings are at least that component's minimum.
Cross-component crossings are nonnegative. Conversely, concatenate optimal
component layouts in one common component order on every species row. No
cross-component link pairs cross. Therefore C* = sum_k C*_k, provided all hard
couplings are included and component placement has no additional constraint.
This gives an unconditional reason to remove a solved component from search.

The Z chain would be a **strand component** (the multispecies analogue of the
Strand/Bundle template) if it has one chromosome per occupied row and is truly
independent. Let its orientation on row i be s_i in {-1,+1}, and let
c_i(s_i,s_(i+1)) be its exact crossing count at that interface. Then:

    D_(i+1)(t) = min_s [D_i(s) + c_i(s,t)].

The recurrence retains both boundary signs, rather than fixing a locally good
choice too soon. It examines at most four transitions per adjacent row pair,
then reconstructs all signs. A 31-row chain need not enumerate 2^31 sign vectors.
The transition costs still require reading/counting all links. The result may
be a positive minimum because within-chromosome rearrangements remain immutable.

### Lossless cell aggregation

A cell contains links with the same two endpoint chromosomes. Every unordered
link pair belongs either within one cell or between two cells:

- Within one cell: compute four sign-dependent inversion counts by sorting and
  a Fenwick counter, retaining strict endpoint ties.
- Two cells sharing one endpoint chromosome: count less-than and greater-than
  coordinate pairs using sorted endpoints and binary search for each sign.
- Two cells sharing neither chromosome: the crossing cost depends only on two
  chromosome precedence bits; its multiplicity is the product of cell sizes.

These cases partition the original link pairs exactly. Summing their tables
reconstructs the same binary objective. This compresses *cost computation*,
not biological evidence. It cannot make the unresolved MILP intrinsically easy,
and a graph with many distinct cells can still be expensive. Timing below does
not claim a real-Pieridae speedup.

The deadline is cooperative, not a hard wall-clock guarantee: parsing, sorting,
presolve, lower-bound preparation, one cell calculation and solver API matrix
construction can overrun a nominal deadline. Slurm remains the hard limit; the
runner checkpoints the presolved layout before entering general construction.

## Phylogeny and the two chromosome organizations

Do not contract all members of a genus or a visually similar karyotype into a
single decision variable without a proof. Similar layouts can still differ at
one consequential endpoint, and optimal pairwise states can conflict at a
shared middle row. Even an exactly solved pair cannot generally be frozen
inside a three- or four-species problem.

The appropriate next mechanism is **separator-state dynamic programming**:
summarize an easy subproblem by a cost table conditional on its boundary states,
not one preferred solution. Combine those tables by min-plus operations. The
exponential parameter is the unresolved boundary state space, not merely the
number of species. Repeated rows can be contracted only when the reduced factor
tables and reconstruction mapping establish equivalence for every boundary
state. Otherwise use coordinated clade-wide proposals and accept them only if
the full unchanged objective improves. Those are proposals, not proof pruning.

For species-row order, use a fixed, supported phylogeny and allow rotations of
children only. Optimal leaf ordering (Bar-Joseph, Gifford & Jaakkola 2001,
doi:10.1093/bioinformatics/17.suppl_1.S22) supplies the relevant endpoint-state
DP for a *fixed pairwise dissimilarity matrix*. It is not automatically a joint
chromosome-layout optimizer. Polytomies must stay explicit; the supplied genus
scaffold records approximate tips and an Anthocharis/Euchloe coarsening, so it
is not a fully verified species topology. Do not silently resolve it to improve
the drawing.

A separate experiment should compare the current order with supported tree
rotations, rebuilding adjacent-pair links from the same complete BUSCO hit
inventory. Existing adjacent-only tables cannot supply new interfaces. Freeze
the universe of underlying evidence, record each displayed pair and its link
count, and report per-interface plus fixed-pair evaluation metrics. A smaller
sum after changing adjacencies is not a certified improvement of the original
fixed-row problem. Inferring where the ancestral/Pieris-like change occurred
would require character mapping/ancestral reconstruction on the tree, beyond
this visualization optimizer.

Relevant work on exact layered crossing minimization and speedups:
Wilson et al. (2024), doi:10.1109/TVCG.2024.3456349,
https://vis.khoury.northeastern.edu/pubs/Wilson2024EvaluatingExtendingSpeedup/.
Its reductions motivate investigation; they require adaptation/proof for
fixed within-chromosome coordinates and whole-chromosome flips here.

## Pegasus test

Use `validation/real_data/submit_pieridae_components.sh EXISTING_VIEW 600`.
For diagnosis without general model construction, append `--audit-only`.
The intended first source is the existing 31-species view:

    <existing-relatedness-run>/topology_views/without-leptidea.SW5Gwr

Its presence and receipts are checked on Pegasus, not assumed by this review.
The wrapper reuses its existing environment and freezes this branch's Python
source, the original view's plotting code when available, and hashes. It uses
one CPU and 16 GB, with a Slurm allowance exceeding the requested search budget.
It never overwrites the previous result. It tests only the canonical all-link
objective and starts from the better full-score layout of the existing pair.

Outputs:

- `audit.json`: component membership/status, full score, exact strand transition
  tables, explicit-Z-label component IDs, old pair-preparation size estimates.
- `chromosome_components.tsv`, `chromosome_pair_support.tsv`: identify every
  Z-component member and all chromosome partners; `chr1` is not assumed to be Z.
- `presolved_layout.tsv`: recoverable output written before general search.
- `optimization/checkpoint_layout.tsv`, `solver.json`, `COMPLETE.json`:
  incumbent, bounds, method and final checks.
- `pdfs/`, if the source view provides its frozen plotter and its required tools.
  This rendering path still needs validation on Pegasus.

First acceptance gate: same row order, input hashes, chromosome set, links and
coordinates; nonincreasing full score; consistent component order; independently
rescored final layout. Second: are the explicit Z component(s) actually strands?
If not, inspect partner counts before any biological filtering decision. Third:
compare preparation time and bound progress at identical budgets. A successful
small regression suite is not a claim that the 31-species optimum will close.

## Verification performed locally

Twelve focused tests (including the runner integration check): exhaustive flips for random 2–5-row strands; a 31-row
strand excluded from MILP construction; adjacent hard equations and nonlocal
deferral; missing-row adjacency; a one-link bridge preventing false separation;
new versus old factor model versus canonical scores over all legal states of
four tiny fixtures; tie/reversal arithmetic; mixed-component solving through
SciPy; honest interruption; legacy/opt-in agreement. Public CLI quickstart:
10 to 0 crossings. Probe audit and solve stages, Bash syntax, mocked submission
and frozen worker execution passed. Python optimization-mode checks also run.
No heavy Pegasus job, real input audit, HighsPy backend regression, full fixture
suite, R installation or real-data rendering was run here.
