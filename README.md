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

The next validation work uses controlled paired presentation tangles and
12/15/20/30-chromosome scaling probes to identify the real scaling wall before
another solver mechanism is added. No solver should be implemented in a way
that contradicts [RULES.md](RULES.md).

## Core formulation

For each species, the legal display state is a constrained signed permutation of whole chromosomes:

```
state = chromosome permutation + whole-chromosome orientation
```

A chromosome may move as a whole or reverse as a whole. Its internal genomic order is immutable.

The intended computational strategy is:

```
homology data
    ↓
multispecies hypergraph / incidence graph
    ↓
connected-component decomposition
    ↓
orientation constraint propagation
    ↓
precedence / interval constraints
    ↓
contract forced structure
    ↓
extract unresolved cyclic kernels
    ↓
exact optimization where feasible
    ↓
spectral initialization only where needed
    ↓
crossing / bundle / displacement minimization
    ↓
re-expand solved structure
    ↓
metrics + visualization
```

## Repository map

- **[RULES.md](RULES.md)** — normative invariants. These define the problem and are protected from casual modification.
- **[RULES_CHANGELOG.md](RULES_CHANGELOG.md)** — explicit record of every approved change to the rule set.
- **[docs/problem_definition.md](docs/problem_definition.md)** — formal statement of the optimization problem.
- **[docs/data_model.md](docs/data_model.md)** — canonical input objects and representations.
- **[docs/pipeline_manifest.md](docs/pipeline_manifest.md)** — stepwise computational plan.
- **[docs/concepts.md](docs/concepts.md)** — mathematical concepts, with status and purpose.
- **[docs/problem_algebra.md](docs/problem_algebra.md)** — emerging signed-permutation/GF(2)/factorized algebraic formulation.
- **[docs/method_stage_audit.md](docs/method_stage_audit.md)** — which graph, DP, and search stages are active, diagnostic, preparatory, or deferred.
- **[validation/benchmark/README.md](validation/benchmark/README.md)** — historical, paired, and scaling benchmarks plus vector reports.
- **[examples/README.md](examples/README.md)** — executable synthetic invariance fixtures.
- **[validation/simulator/README.md](validation/simulator/README.md)** — forward chromosome-evolution simulator, deliberate display-tangle induction, and visual validation.
- **[docs/secondary_ancestral_inference.md](docs/secondary_ancestral_inference.md)** — post-core benchmark plan for ancestral fusion/fission/inversion history inference.

## Design principle

A mathematical concept belongs in SynTangle only if it helps us **represent the data, remove impossible choices, solve the remaining legal choices, or quantify the result**. Interesting mathematics that does not improve one of those jobs remains experimental or deferred.
