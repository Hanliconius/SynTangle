# Algebraic formulation of the SynTangle problem

SynTangle is not yet presented as a finished algebraic theory, but the current
solver already has a natural algebraic structure. Making that structure
explicit is useful because it separates biological invariants from legal
display transformations and clarifies why decomposition and monotone reduction
are valid.

## 1. Legal state space

For species s with n_s chromosomes, a display state consists of:

- a permutation of complete chromosomes;
- an orientation bit for each complete chromosome.

Ignoring additional hard constraints for the moment, the legal transformations
form the signed-permutation (hyperoctahedral) group

B_n = C2^n semidirect-product S_n.

Here S_n reorders whole chromosomes and C2^n flips whole chromosomes. There is
no generator corresponding to a subchromosomal move, which encodes R1--R4
directly into the state space.

For multiple species the unconstrained display space is the direct product of
the species-specific signed-permutation groups.

Observed homology, block coordinates, and within-chromosome order are evidence
used to score states in that space; they are not themselves generators of the
transformation group.

## 2. Orientation constraints as linear algebra over GF(2)

Write each chromosome orientation as a bit x_i in GF(2). A hard relative
orientation relation has the form x_i XOR x_j = b_ij.

The full set of hard orientation equations can therefore be written A x = b
modulo 2. When consistent, the legal orientation solutions are an affine space
x = x0 + ker(A).

OrientationBasis is the computational representation of this object:

- base_assignment is x0;
- free_flip_groups span the remaining null-space degrees of freedom.

Stage 15 deliberately keeps this compact basis rather than expanding it into
all 2^k complete assignments.

## 3. Component factorization

The chromosome--homology incidence graph identifies disconnected evidence
components. Under R7, if no objective term connects two components, the legal
state space factors as a direct product over components and the crossing
objective is additive across those components.

This is why connected-component decomposition is not just a visualization
convenience: it is an exact factorization of the optimization problem.

Further decomposition into bridges, articulation points, biconnected blocks,
2-cores, or structural kernels is valid for optimization only when the
remaining objective also factorizes over the proposed pieces. At present those
decompositions are computed and visualized, but Stage 15 does not yet assume
that every structural kernel is an independent optimization factor.

## 4. Crossing objective as a chain factor graph

For species layers in order 1 through m, the current crossing score is a sum
over adjacent-layer costs:

C(g_1,...,g_m) = sum_s C_(s,s+1)(g_s,g_(s+1)).

For a fixed orientation assignment this is a chain-structured factor graph.
That factorization is what makes the exact species-layer dynamic program
possible.

The one-layer subset dynamic program solves a related conditional problem:
given neighboring layers, choose the chromosome order of one species that
minimizes its incident pairwise crossing contribution.

## 5. Equivalence, quotienting, and normal forms

Several legal states may be representationally equivalent. Disconnected
components, for example, can be placed in one common canonical order without
changing the biological interpretation.

This suggests viewing canonicalization and contraction as quotient operations:
states related by a proven equivalence relation are represented by one
canonical element of the quotient space.

The constrained optimum C* defines a set of optimal layouts rather than
necessarily one unique layout. A deterministic tie-break can choose one
representative as a display normal form, while the full equivalence/optimal
class remains scientifically important.

## 6. Monotone residual search

R17 gives the search a monotone algebraic interpretation. At each stage the
solver replaces the current feasible domain with a subset or quotient that
contains every state still capable of improving the incumbent:

D0 contains D1 contains D2 contains ...

A variable may disappear because it is:

- forced by hard equations;
- equivalent to an already represented choice;
- excluded by a valid lower bound;
- irrelevant once lower and upper bounds meet.

The original biological evidence is not quotiented away. Only the unresolved
display decision space shrinks.

## 7. What would make this a fuller algebra of chromosome untangling?

A more complete theory would explicitly define:

1. the generators and relations of legal whole-chromosome transformations;
2. equivalence relations induced by disconnected or symmetric structure;
3. composition rules for component solutions;
4. a canonical normal-form convention for non-unique optima;
5. admissible contractions and their inverses;
6. lower-bound-preserving maps from the full problem to relaxed problems;
7. eventually, tree-constrained evolutionary-event operators for the separate
   ancestral-inference problem.

It is therefore accurate to say that SynTangle is developing an algebraic
description of the problem, but not yet that it has a finished named algebra.
The strongest existing pieces are signed permutations, GF(2) orientation
linear algebra, direct-product component factorization, chain-factor objective
decomposition, quotient/canonical-state reasoning, and monotone restriction of
the residual feasible domain.
