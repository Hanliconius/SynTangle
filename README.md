# SynTangle

**SynTangle** is a framework for producing minimally tangled multispecies chromosome-synteny layouts without changing the biological order encoded by the input assemblies.

> Intrinsic tangledness is biology. Layout tangledness is noise.

SynTangle separates those two things. It treats whole chromosomes as the movable units, preserves true within-chromosome order exactly, decomposes the multispecies correspondence graph as far as possible, and optimizes only the unresolved legal ordering/orientation choices that remain.

## Current status

Implementation has begun with **executable validation infrastructure**: machine-readable ground-truth fixtures plus a recovered/cleaned forward chromosome simulator for fusion, fission and inversion, deliberate whole-chromosome tangle induction, and visual audit output.

The SynTangle solver itself is still being built from the specification outward. No solver should be implemented in a way that contradicts [RULES.md](RULES.md).

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
- **[examples/README.md](examples/README.md)** — executable synthetic invariance fixtures.
- **[validation/simulator/README.md](validation/simulator/README.md)** — forward chromosome-evolution simulator, deliberate display-tangle induction, and visual validation.
- **[docs/secondary_ancestral_inference.md](docs/secondary_ancestral_inference.md)** — post-core benchmark plan for ancestral fusion/fission/inversion history inference.

## Design principle

A mathematical concept belongs in SynTangle only if it helps us **represent the data, remove impossible choices, solve the remaining legal choices, or quantify the result**. Interesting mathematics that does not improve one of those jobs remains experimental or deferred.
