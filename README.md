# SynTangle

**SynTangle** is a framework for producing minimally tangled multispecies chromosome-synteny layouts without changing the biological order encoded by the input assemblies.

> Intrinsic tangledness is biology. Layout tangledness is noise.

SynTangle separates those two things. It treats whole chromosomes as the movable units, preserves true within-chromosome order exactly, decomposes the multispecies correspondence graph as far as possible, and optimizes only the unresolved legal ordering/orientation choices that remain.

## Current status

SynTangle now has an executable exact/bounded solver stack: canonical fixture
parsing, chromosome↔homology incidence components, GF(2) orientation
propagation, exact species-layer dynamic programming, exact one-layer subset
DP, monotone residual branch-and-bound, seeded local-search incumbents, and
auditable initial-vs-optimized visualization. The original 27-case simulation
benchmark is solved to proven optimum throughout.

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

Stage 18 also introduces process-level parallelism across exact independent
incidence components in monotone branch-and-bound and an explicit residual
variable/factor graph for measuring the finer decision structure that remains
after legal reductions. Finer residual-factor/separator parallelism is the next
solver integration step; it will be activated only when objective independence
is proved. No solver should contradict [RULES.md](RULES.md).

## Core formulation

For each species, the legal display state is a constrained signed permutation of whole chromosomes:

```
state = chromosome permutation + whole-chromosome orientation
```

A chromosome may move as a whole or reverse as a whole. Its internal genomic order is immutable.

The current computational strategy is:

```
multispecies homology hypergraph
    ↓ exact sparse incidence representation
independent incidence components
    ↓
GF(2) orientation propagation + canonicalization
    ↓
residual order/orientation decision-factor graph
    ↓
exact layer DP where feasible
    ↓ otherwise
monotone branch-and-bound with legal incumbents/bounds
    ↓
proven optimum OR explicit lower/upper-bound gap
    ↓
visual audit + complexity diagnostics
```

Independent incidence components may be solved in separate worker processes.
Disconnected residual-factor pieces and small separators are the next
proof-preserving decomposition layer.

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
- **[validation/benchmark/README.md](validation/benchmark/README.md)** — historical, paired, and scaling benchmarks plus vector reports.
- **[examples/README.md](examples/README.md)** — executable synthetic invariance fixtures.
- **[validation/simulator/README.md](validation/simulator/README.md)** — forward chromosome-evolution simulator, deliberate display-tangle induction, and visual validation.
- **[docs/secondary_ancestral_inference.md](docs/secondary_ancestral_inference.md)** — post-core benchmark plan for ancestral fusion/fission/inversion history inference.

## Design principle

A mathematical concept belongs in SynTangle only if it helps us **represent the data, remove impossible choices, solve the remaining legal choices, or quantify the result**. Interesting mathematics that does not improve one of those jobs remains experimental or deferred.
