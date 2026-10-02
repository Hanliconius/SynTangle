# Stage 0 executable fixtures

These JSON files are the first machine-readable specification of SynTangle's behavior.

They are deliberately **synthetic**. Every fixture has a known answer by construction, so a failed assertion can be attributed to the implementation rather than to uncertain biology.

## Representation

Each fixture contains:

- species;
- whole chromosomes;
- immutable, coordinate-ordered homologous block occurrences;
- optional explicit whole-chromosome orientation constraints;
- an `expected` object containing facts that future code must recover.

The fixture format is described by `fixture.schema.json`.

The chromosome `display_rank` is deliberately distinct from genomic block order. Display rank may change. Block coordinate order may not.

## Current fixtures

1. **perfect_1to1_30x3.json** — 30 independent 1:1:1 components; no hard kernel.
2. **fusion_chain_tree.json** — connected fusion component whose incidence graph is still a tree.
3. **fusion_chain_closed_cycle.json** — closes one independent cycle and creates a cyclic kernel.
4. **balanced_orientation_cycle.json** — orientation equations are solvable by propagation.
5. **frustrated_orientation_cycle.json** — one irreducible signed inconsistency.
6. **forbidden_subchromosomal_pretty_solution.json** — zero crossings would require an illegal internal reorder.
7. **layout_noise_only.json** — three crossings are purely display-induced and can be removed by whole-chromosome ordering.

## Rule for future fixtures

A fixture should be small enough that its expected result can be understood by inspection. Large or realistic datasets belong in the real-data validation tier, not here.

These files should eventually be run in CI and remain fast enough for every pull request.
