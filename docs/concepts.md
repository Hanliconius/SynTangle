# Mathematical concepts and their jobs

A concept belongs in SynTangle only if it has a clear computational or descriptive job.

## CORE

### Multispecies hypergraph

**Job:** faithfully represent one homology relationship shared across more than two species without all-pairs duplication.

**Why retained:** pairwise expansion loses the natural multiway object and can scale quadratically with species occupancy.

### Incidence graph

**Job:** expose the hypergraph to standard graph algorithms.

**Used for:** connected components, bridges, articulation points, 2-core, biconnected components, cycle basis, cycle rank, and kernel decomposition.

### Signed graph / binary algebra over `GF(2)`

**Job:** solve whole-chromosome orientation.

Each relationship gives (x_i\oplus x_j=b_{ij}). Balanced components are solved by propagation; inconsistent cycles identify frustration.

### Constrained signed permutations

**Job:** define the actual legal layout state.

Whole chromosomes may permute and reverse. Chromosome interiors do not.

The unconstrained algebraic state resembles a signed permutation group, but the biological data restrict the search to an admissible subset (\Omega).

### Crossing / inversion minimization

**Job:** quantify and minimize layout tangledness.

For fixed ordered layers, crossings correspond closely to inversions. The constrained minimum is the central visual objective.

## ADOPTED / HIGH PRIORITY

### Graph decomposition

**Job:** make the hard problem small before solving it.

Connected components, bridges, articulation points, 2-cores, biconnected components, and contraction of forced structure are algorithmic machinery, not merely descriptive metrics.

### Cycle rank and cycle basis

**Job:** describe independent cyclic constraint and locate the parts of the graph where ambiguity/conflict can persist.

SynTangle keeps two related views distinct:

- the **raw incidence graph**, which retains every homologous block/anchor and therefore reflects evidence density as well as structure;
- a **chromosome-signature structural projection**, which bundles homology groups having the same chromosome-incidence multiset before computing chromosome-structural cycle/core metrics.

The structural projection is only a decomposition aid. It never replaces the full ordered block data used for crossings, inversions, duplications, or within-chromosome order. Copy multiplicity within a homology group remains explicit.

Do not enumerate all simple cycles in large graphs.

### Precedence / consecutive-order constraints

**Job:** represent order that is already forced by immutable chromosome interiors.

### PQ-trees / PC-trees or related machinery

**Job:** compactly represent families of legal orders satisfying consecutive/reversible constraints.

**Status:** high-priority to test. Particularly attractive because a chromosome behaves like a reversible ordered chain: native order or complete reversal, never internal permutation.

### Frustration index

**Job:** quantify irreducible inconsistency among whole-chromosome orientation constraints.

**Status:** retain as a metric and as a way to isolate a small conflicting orientation kernel. Exact optimization may be replaced by bounds/heuristics for large cases.

## EXPERIMENTAL

### Spectral graph / hypergraph ordering

**Job:** provide a good continuous initialization when the unresolved legal kernel remains too large for exact search.

**Important restriction:** spectral coordinates may guide placement of whole chromosomes/supernodes only. They may never reorder internal chromosome content.

### Bundle/block crossing metrics

**Job:** capture perceptual tangledness when many pairwise crossings visually form one bundled crossing event.

### Multiscale tangledness

**Job:** ask at what genomic resolution complexity appears while always preserving true within-chromosome order.

## DEFERRED

### Braid theory

**Possible job:** describe residual interweaving after a legal optimal/near-optimal ordering is found.

**Why demoted:** SynTangle primarily cares about final constrained chromosome order, not the path by which swaps were performed. Signed permutations/inversion length capture the central ordering algebra more directly. Braid analysis can return later if it adds biological information not already captured by crossings/permutations.

## NOT PART OF THE CORE MODEL

- unrestricted subchromosomal spectral seriation;
- arbitrary block flips or swaps;
- full braid-group optimization as the primary solver;
- one undifferentiated “complexity score”;
- eager all-pairs expansion of every multispecies hyperedge.

These may be mathematically interesting but either violate the biological invariants or do not currently buy enough computational value.
