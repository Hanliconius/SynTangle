# Formal problem definition

## 1. Biological object

For species (s), let its chromosomes be

[
mathcal C_s = \{c_{s1}, c_{s2}, \ldots, c_{sn_s}\}.
]

Each chromosome is an immutable ordered interval containing observed homologous block occurrences. An occurrence can be represented as

[
b=(s,c,start,end,strand,h,confidence),
]

where (h) is a multispecies homology-group identifier.

Internal order along a chromosome is data. It is never an optimization variable.

## 2. Legal state of one species

A display state for species (s) is

[
g_s=(\pi_s,\epsilon_s),
]

where:

- (\pi_s\in S_{n_s}) is a permutation of **whole chromosomes**;
- (\epsilon_s\in\{-1,+1\}^{n_s}) is the orientation of each whole chromosome.

For a chromosome of length (L), a coordinate (x) transforms as

[
x'=
\begin{cases}
x,&\epsilon=+1\\
L-x,&\epsilon=-1.
\end{cases}
]

No other internal transformation is legal.

Across all species, the state is

[
G=(g_1,g_2,\ldots,g_S).
]

The biological constraints define an admissible state set

[
\Omega \subseteq \prod_s (C_2^{n_s}\rtimes S_{n_s}).
]

SynTangle optimizes only over (\Omega), never over unrestricted subchromosomal permutations.

## 3. Homology representation

A homology group connecting occurrences in multiple species is naturally represented as a hyperedge. The complete dataset is therefore a multispecies hypergraph.

For graph decomposition and most algorithms, SynTangle may use the corresponding bipartite incidence graph:

```
chromosome/block occurrence ── homology group ── chromosome/block occurrence
```

The hypergraph is the faithful biological representation; the incidence graph is an algorithmic projection.

## 4. Orientation constraints

Whole-chromosome orientation can be represented by binary variables

[
x_i\in\{0,1\}.
]

A homologous orientation relationship between chromosomes (i) and (j) imposes

[
x_i\oplus x_j=b_{ij},
]

with (b_{ij}=0) for concordant and (b_{ij}=1) for reversed orientation.

A signed cycle is consistent when its XOR sum is zero. In a balanced component, all orientations are determined up to global reversal once one chromosome state is fixed.

## 5. Ordering constraints

Observed chromosome interiors impose immutable order/precedence constraints on the homologous material they contain. A chromosome may be reversed as a whole, but its internal sequence may not be re-permuted.

The legal permutation space should therefore be represented as compactly as possible using precedence, interval/consecutive, component, and reversible-chain constraints rather than by enumerating all (n!\) orders.

## 6. Objective

For a legal state (G\in\Omega), define one or more display costs.

Primary:

[
C(G)=\text{weighted number of synteny-link crossings}.
]

Additional candidates:

[
B(G)=\text{bundle/block crossing cost},
]

[
D(G)=\text{weighted positional displacement from a monotone/diagonal arrangement}.
]

Version 1 should keep these quantities separate or combine them only through an explicit documented lexicographic/weighted objective.

The constrained optimum is

[
C^*=\min_{G\in\Omega}C(G).
]

For a drawing (G),

[
T_{excess}(G)=C(G)-C^*.
]

Thus (T_{excess}=0) means the drawing is minimally tangled under the biological constraints, while (C^*) measures residual/intrinsic crossing complexity.

When exact global optimality is computationally unavailable, SynTangle must distinguish **best known** from **proven optimum**.

## 7. Computational philosophy

Do not solve a large permutation problem when graph structure proves that most choices are irrelevant or forced.

The intended reduction is:

```
full input
  → connected components
  → orientation propagation
  → precedence/interval constraints
  → contract forced structure
  → bridge/articulation/core decomposition
  → unresolved kernel
  → optimize kernel only
  → expand
```

Practical difficulty is expected to depend more strongly on the size/structure of the largest unresolved kernel than on the raw number of chromosomes.
