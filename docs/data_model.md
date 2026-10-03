# Canonical data model

The initial implementation should normalize upstream synteny sources into a small set of explicit tables/objects. Upstream file formats are adapters; the internal representation should not depend on GENESPACE, MCScanX, a particular aligner, or a particular plotting package.

## Chromosome

Required fields:

```text
species_id
chromosome_id
length
input_order            # optional original display order; not biological order within chromosome
input_orientation      # optional +1 / -1 original display orientation; not biological
metadata               # optional
```

## Homologous block occurrence

One physical occurrence of one homology group:

```text
homology_id
species_id
chromosome_id
start
end
strand                 # relative orientation where defined
confidence             # optional
source                  # optional provenance
occurrence_id           # unique
```

The ordering of occurrences is derived from genomic coordinates and is immutable.

## Homology group

A multispecies relationship:

```text
homology_id
occurrence_ids[]
weight                  # default may derive from homologous length
confidence              # optional group-level confidence
metadata
```

A homology group is represented as one hyperedge, not eagerly expanded into every species pair.

## Derived chromosome correspondence

For species (a,b), chromosome correspondence may be summarized on demand as

[
S^{ab}_{ij}=\sum_{h\in c_i^a\leftrightarrow c_j^b} w_h.
]

These pairwise matrices are projections/caches, not the primary data representation.

## Layout state

For each species:

```text
chromosome_order[]
chromosome_orientation{}   # +1 / -1
```

There is no field for subchromosomal reordering because that operation is illegal.

## Constraint state

Derived constraints should be first-class and inspectable:

```text
connected_component_id
orientation_equations
precedence_constraints
contiguity / reversible-chain constraints
contracted_supernodes
bridge / articulation decomposition
2-core membership
biconnected_component_id
cycle_basis
unresolved_kernel_id
```

## Audit record

Every solved layout should retain:

```text
input_hash / provenance
initial_state
final_state
whole_chromosome_moves
whole_chromosome_flips
constraints_used
solver / strategy
random_seed, if any
objective_before
objective_after
optimality_status       # proven optimum / bounded / best known
metrics
```
