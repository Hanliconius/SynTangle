# Secondary goal: ancestral structural-event inference

This is a **post-core research objective**, not part of the current SynTangle untangling problem.

Once SynTangle can recover a minimally tangled legal representation of extant chromosome homology, a modified version of the same graph/decomposition machinery may be useful for asking a different question:

> Which ancestral sequence(s) of chromosome fusion, fission and inversion events could plausibly have generated the observed extant chromosome structures?

## Keep the two problems separate

The core SynTangle problem is descriptive/layout-oriented:

```text
extant chromosome homology
        ↓
legal whole-chromosome order/orientation
        ↓
minimally tangled extant representation
```

The proposed secondary problem is historical:

```text
extant chromosome homology
        ↓
candidate ancestral chromosome states
        ↓
candidate fusion/fission/inversion histories
        ↓
ranked or bounded set of compatible histories
```

A visually simple extant layout is not itself evidence that one particular evolutionary history occurred. The ancestral-inference layer must therefore remain separate from the untangling objective and must explicitly represent non-identifiability.

## Why the current architecture may help

Several objects already needed for SynTangle are also useful for structural-history inference:

- multispecies homology groups;
- chromosome↔homology incidence components;
- immutable within-chromosome block order;
- orientation constraints;
- bridges, articulation points and 2-cores;
- cycle rank and a cycle basis;
- fusion/fission-like multi-chromosome correspondence patterns.

The decomposition may let ancestral inference operate component-by-component rather than searching an unconstrained global event history.

## First benchmark: hidden-truth simulation

The validation simulator is an unusually clean benchmark because it records the true evolutionary event sequence separately from the presentation-only tangle log.

For an ancestral-inference benchmark:

1. simulate an ancestor;
2. generate descendants with logged fusion/fission/inversion events;
3. discard the evolution log from the inference input;
4. optionally tangle only the display state;
5. run SynTangle to obtain the legal extant representation;
6. infer one or more candidate ancestral event histories from extant data alone;
7. reveal the hidden evolution log only for scoring.

The presentation-only `tangle_log` must never count as evolutionary truth.

## What success should mean

Do not score only exact event-string recovery. Distinct event sequences can generate the same extant chromosome state.

Useful benchmark quantities include:

- recovery of ancestral chromosome adjacencies;
- recovery of event types and counts;
- recovery of affected chromosome/homology components;
- breakpoint error/tolerance for inversions and fissions;
- edit distance between inferred and simulated event sequences where order is identifiable;
- whether the true history lies in the returned compatible/near-optimal set;
- size of the equivalence class of equally supported histories;
- confidence or support for identifiable events;
- runtime as structural complexity and species number increase.

## Candidate inference formulations for later testing

These are research directions, not commitments:

- minimum-event/parsimony search on chromosome states;
- dynamic programming on decomposed components;
- branch-and-bound over structural-event histories;
- likelihood/Bayesian scoring if an explicit event-rate model is justified;
- ancestral adjacency reconstruction followed by event decomposition;
- reuse of cycle/core structure to localize historically ambiguous regions.

A single maximum-parsimony history should not be reported as uniquely true when multiple histories are observationally equivalent.

## Development gate

Do not start this module until the extant untangling pipeline can:

- preserve all SynTangle invariants;
- reproduce the synthetic fixture expectations;
- separate intrinsic from layout tangledness;
- produce stable chromosome order/orientation output;
- expose its decomposition and audit record.

At that point the simulator can support a second benchmark suite specifically for historical reconstruction.
