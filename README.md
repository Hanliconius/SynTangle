# SynTangle

**SynTangle** is a framework for producing minimally tangled multispecies chromosome-synteny layouts without changing the biological order encoded by the input assemblies.

> Intrinsic tangledness is biology. Layout tangledness is noise.

SynTangle separates those two things. It treats whole chromosomes as the movable units, preserves true within-chromosome order exactly, decomposes the multispecies correspondence graph as far as possible, and optimizes only the unresolved legal ordering/orientation choices that remain.

## Current status

The Pegasus development branch also contains experimental joint MILP and hybrid
neighborhood solvers. See the [benchmark history](docs/benchmark_history.md) for
methods, tested commits, job IDs, results and limitations, and the
[archived global-method comparison](validation/benchmark/results/2026-10-06-global-methods/summary.md)
for the complete collector output. The latest hybrid comparison reports matching global bounds on all nine stress
presentations, including the largest cases at C=8377. These experiments have not
replaced the default production solver; raw proof/audit review remains outstanding.

SynTangle now has an executable exact/bounded solver stack: canonical fixture
parsing, chromosome↔homology incidence components, GF(2) orientation
propagation, recursive exact residual factor elimination, exact species-layer
dynamic programming as an independent fallback, exact one-layer subset DP,
monotone branch-and-bound, seeded local-search incumbents, and auditable
initial-vs-optimized visualization. The original 27-case simulation benchmark
is solved to proven optimum throughout.

Structural projection, bridges/articulation points, 2-cores, biconnected
blocks, and the fundamental cycle basis are still computed and visualized, but
they are currently diagnostic rather than being allowed to prune the optimizer
without a proof that the remaining crossing objective factorizes accordingly.

Current validation now includes controlled paired presentation tangles,
historical chromosome-scale probes, a coupled high-complexity stress ladder,
and an orthogonal biologically grounded benchmark that varies species count,
chromosome count, and rearrangement burden separately. The orthogonal profile
uses independent extant descendants from one hidden ancestor so increasing
species count does not also increase cumulative rearrangement depth. Long
benchmark runs checkpoint every completed case and can resume.

Stage 18 introduced the explicit residual variable/factor graph and
process-level parallelism across independent incidence components. Stage 20
makes that residual graph active solver machinery: exact factor tables contract
variables that prove irrelevant, primal-graph leaves are eliminated,
disconnected residual pieces are solved separately, small articulation
variables are conditioned, and min-fill elimination is used only for the
irreducible remainder. No solver should contradict [RULES.md](RULES.md).

## Current development direction

The latest Pegasus experiments shift the emphasis from refining branch-search
bounds to optimizing chromosome order and orientation jointly across species.
Joint MILP reported matching layout scores and global bounds at C=341 on the
6-species stress cases and C=1258 on the 8-species cases. On the 10-species
cases, multi-species neighborhood search found better layouts than the tested
MILP and pipeline runs, but the global proof gaps remain large. See the
[benchmark history](docs/benchmark_history.md) for matched results and caveats.

The next candidate combines these strengths:

- Use incidence components, GF(2) orientation propagation and crossing factors
  to construct a joint order/orientation model.
- Improve a feasible layout with bounded multi-species neighborhoods, then pass
  it as a solver start to global MILP.
- Reuse sparse model matrices, retain proven-zero components, and allocate
  unused time to unresolved components rather than restart unchanged work.
- Preserve every applicable proven reduction across stages. A heuristic layout
  supplies an upper bound; its preferred ordering does not by itself justify
  excluding alternatives. Neighborhood restrictions are conditional and are
  released before global search.

This remains a graph-informed approach. The joint formulation addresses
ordering choices shared by neighboring species that independent conditional
optimizations can miss. Structural diagnostics still require a proof before
they can remove feasible choices.

The [hybrid comparison](validation/benchmark/results/2026-10-06-hybrid-methods/summary.md)
now reports all nine stress presentations solved to matching global bounds at
600-second allowances. The three largest cases reach C=8377 in 428–501 seconds;
adaptive neighborhoods yield better layouts at shorter budgets, while plain
MILP also reaches the final optimum. Stage-by-stage efficiency attribution and
the setup strict-improvement audit output still require review. Graph/GF(2) reductions and model
data are reused in the hybrid, but importing every residual-domain reduction
and reconstruction mapping from the older exact backends has not yet been
established. That handoff contract remains a requirement before claiming full
reduction continuity or promoting the experimental backend.

## Core formulation

For each species, the legal display state is a constrained signed permutation of whole chromosomes:

```
state = chromosome permutation + whole-chromosome orientation
```

A chromosome may move as a whole or reverse as a whole. Its internal genomic order is immutable.

The default solver's computational strategy is:

```
multispecies homology hypergraph
    ↓ exact sparse incidence representation
independent incidence components
    ↓
GF(2) orientation propagation + canonicalization
    ↓
residual order/orientation decision-factor graph
    ↓
contract → split → leaf-eliminate → condition separators → repeat
    ↓
exact min-fill elimination on surviving kernels
    ↓ if configured residual caps are exceeded
exact layer DP
    ↓ otherwise
monotone branch-and-bound with legal incumbents/bounds
    ↓
proven optimum OR explicit lower/upper-bound gap
    ↓
visual audit + complexity diagnostics
```

Independent incidence components may be solved in separate worker processes.
Disconnected residual-factor pieces and one-variable articulation separators
are now active exact decomposition layers. Larger multi-variable separators are
deferred until a surviving benchmark kernel demonstrates that they are needed.

## Repository map

- **[RULES.md](RULES.md)** — normative invariants. These define the problem and are protected from casual modification.
- **[RULES_CHANGELOG.md](RULES_CHANGELOG.md)** — explicit record of every approved change to the rule set.
- **[docs/problem_definition.md](docs/problem_definition.md)** — formal statement of the optimization problem.
- **[docs/data_model.md](docs/data_model.md)** — canonical input objects and representations.
- **[docs/pipeline_manifest.md](docs/pipeline_manifest.md)** — stepwise computational plan.
- **[docs/concepts.md](docs/concepts.md)** — mathematical concepts, with status and purpose.
- **[docs/problem_algebra.md](docs/problem_algebra.md)** — emerging signed-permutation/GF(2)/factorized algebraic formulation.
- **[docs/method_stage_audit.md](docs/method_stage_audit.md)** — which graph, DP, and search stages are active, diagnostic, preparatory, or deferred.
- **[docs/verification_roadmap.tex](docs/verification_roadmap.tex)** — concise narrative of the minimum benchmark/validation chain needed to road-test the method.
- **[docs/benchmark_history.md](docs/benchmark_history.md)** — Pegasus experiment chronology, evidence, limitations and current decisions.
- **[validation/benchmark/README.md](validation/benchmark/README.md)** — historical, paired, and scaling benchmarks plus vector reports.
- **[examples/README.md](examples/README.md)** — executable synthetic invariance fixtures.
- **[validation/simulator/README.md](validation/simulator/README.md)** — forward chromosome-evolution simulator, deliberate display-tangle induction, and visual validation.
- **[docs/secondary_ancestral_inference.md](docs/secondary_ancestral_inference.md)** — post-core benchmark plan for ancestral fusion/fission/inversion history inference.

## Design principle

A mathematical concept belongs in SynTangle only if it helps us **represent the data, remove impossible choices, solve the remaining legal choices, or quantify the result**. Interesting mathematics that does not improve one of those jobs remains experimental or deferred.
