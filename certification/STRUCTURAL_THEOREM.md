# Structural certificates for whole-chromosome synteny layouts

**Status:** proved elementary statements under the precise model below;
experimental implementation with exhaustive small-instance tests. Not independently
peer-reviewed, not a claim of a new theorem in signed graph theory, and not a
claim that the bounds will be tight on our real datasets.

This note adapts standard signed-graph balance and cycle-obstruction ideas to
SynTangle's exact unweighted strict midpoint crossing objective. The relevant
graph is a pairwise projection of the richer ordered homology/factor structure.
We do not assert that an unordered incidence hypergraph determines the optimum.

## Model

Species rows are fixed. A legal layout permutes whole chromosomes within each
row and gives each chromosome a bit f_a, where one means a whole-chromosome flip.
Internal coordinates and observed homologies are immutable. Hard equations are
f_a XOR f_b = p. One-to-one observed homologies join adjacent rows. Count a pair
of links as crossing exactly when their strict endpoint orders disagree; pairs
with a tie on either side count zero.

For each chromosome pair e=(a,b) on adjacent rows, collect all links between a
and b into a **bundle**. Among pairs of links inside that bundle, let:

* P_e = the number with concordant native endpoint orders;
* D_e = the number with discordant native endpoint orders.

Tied pairs contribute to neither count. Counts use exact rational input coordinates.
Each eligible link pair belongs to one bundle, so nothing is counted twice.

## Lemma 1: exact decomposition

For every legal layout l,

    C(l) = sum_e c_e(f_a XOR f_b) + R(l),
    c_e(0)=D_e,   c_e(1)=P_e,   R(l)>=0.

Here R counts all crossing link pairs outside a single bundle. These include
pairs sharing one endpoint chromosome, and pairs sharing neither.

**Proof.** Partition all link pairs into those lying in one bundle and all
remaining pairs. A bundle's endpoint comparisons depend only on internal
coordinates and the two flips, never chromosome positions. If flip parity is
zero, both comparisons retain their native agreement/disagreement; if one end
flips, every strict comparison reverses and concordance becomes discordance.
Each remaining link pair contributes either zero or one. Summing proves the identity.

## Theorem 1: a universal structural lower bound and attainment certificate

Define B0 = sum_e min(P_e,D_e). For every legal layout C(l)>=B0.
If a validated candidate has C(l)=B0, then it is globally optimal.

**Proof.** Each bundle contribution is at least its smaller possible value and
R is nonnegative. A candidate attaining the lower bound attains the global minimum.
This argument uses neither exhaustive search nor a numerical optimization solver.

The count is useful even for positive optima. Example: on one chromosome per
species, three homologies occur as a,b,c in one row and a,c,b in the other.
Native discordance is one and concordance is two. Whole-chromosome flips cannot
improve on one crossing; a native layout attains it, so C*=1.

## Theorem 2: signed cycle penalties strengthen the bound

For each bundle with P_e != D_e, define the preferred parity q_e:
q_e=0 when D_e<P_e, and q_e=1 when P_e<D_e. Set w_e=abs(P_e-D_e)>0.
Then

    c_e(f_a XOR f_b)
      = min(P_e,D_e) + w_e * [f_a XOR f_b != q_e].

Create a signed multigraph whose vertices are chromosomes and whose soft edges
are these bundle preferences. Add hard orientation equations as hard parity
edges. A cycle whose edge parities XOR to one is inconsistent: it is impossible
to satisfy all its equations simultaneously.

For any collection of inconsistent cycles K such that **no soft edge appears
in more than one selected cycle**, define

    B_K = B0 + sum_{Q in K} min_{soft e in Q} w_e.

Then C(l)>=B_K for every legal layout. Hard edges may be shared between cycles.
A cycle containing only hard edges implies that no legal layout exists.
A validated candidate with C(l)=B_K is globally optimal.

**Proof.** XORing the actual endpoint parities around a cycle yields zero,
because each vertex bit appears twice. In a cycle whose prescribed parities XOR
to one, at least one edge must therefore violate its prescribed parity. Legal
layouts cannot violate hard edges, so at least one soft edge pays w_e, at least
the minimum soft weight on the cycle. Since selected cycles share no soft edges,
their mandatory payments cannot double-count the same cost. Add these payments
to B0 and apply Lemma 1. The attainment conclusion follows as in Theorem 1.

Cycle packing is only a sufficient lower bound. The implementation greedily
selects fundamental cycles of a parity forest; it does not find a maximum packing
or generally solve the weighted frustration problem. Different packings can give
different bounds without affecting validity.

## Theorem 3: a constructively solvable structural class

Build the chromosome incidence graph using **all scored links**, ignoring their
orientation costs. Suppose each connected component contains at most one
chromosome from any species row (call this the row-thin condition). Suppose,
further, that all strict bundle preferences together with all hard orientation
equations are consistent.

Then C*=B0, and an optimal layout can be constructed by:

1. Solving the parity equations by propagation.
2. Choosing any total order of incidence components.
3. Ordering the chromosomes in every species by that same component order.

**Proof.** Consistency gives orientations attaining each bundle's minimum.
Consider a pair of links on one adjacent-row interface. If both are in the same
incidence component, the row-thin condition implies they lie in the same bundle.
If they belong to different components, their endpoint chromosomes appear in
the same component order on both rows, so they cannot cross. Consequently R=0
and C=B0. Theorem 1 proves optimality.

Chromosomes lacking scored links are singleton components and contribute no
crossings. Missing rows within a component do not invalidate the proof. Hard
equations may join separate incidence components, provided the combined parity
system remains consistent. Equal-cost bundles impose no orientation preference.

After bundle counts are prepared, consistency and component membership use graph
traversals, without factorial permutation enumeration. The reference implementation
counts bundle pairs quadratically; faster inversion counting can compute P,D later.
This class includes arbitrarily many rows and chromosomes, and permits positive
unavoidable crossings within bundles. It excludes incidence components that
contain multiple chromosomes in a row, a common feature of fusion/fission datasets.

## Important counterexample: balance is not enough

Take two rows with two chromosomes each, and one unique link for every chromosome
pair: the incidence graph is K_2,2. Every bundle contains only one link, so P=D=0,
B0=0, and the preference graph has no inconsistent cycles.

Yet every chromosome ordering has a crossing: the two links from the earlier
chromosome to the later opposite chromosome and from the later chromosome to
the earlier opposite chromosome cross, regardless of flips. In a small explicit
fixture the minimum is exactly one.

Thus a zero orientation-frustration bound does not imply a zero full optimum.
The remaining obstruction is chromosome ordering, represented by R. Neither
ordinary graph acyclicity nor balance alone is claimed to characterize our entire
solvable class. A failed attainment check means 'not certified by this theorem',
not 'unsolvable' and not 'nonoptimal'.

## Research meaning

This separates unavoidable internal-order discordance, unavoidable orientation
frustration, and remaining chromosome-order crossings. It supplies an exact,
solver-independent sufficient certificate for any graph when a layout attains
the bound, plus an efficiently constructible special class.

The next empirical question is how much of the observed residual in planarian
and Leptidea layouts this bound explains. If it is weak, the next theorem target
is an ordering-obstruction bound for R, with a disjointness argument that prevents
double counting. A strong result on a broader structural class would need more
than the unsigned incidence topology.

## Literature grounding

Signed-graph balance and switching are established theory. Zaslavsky's
[Signed Graphs and Geometry](https://people.math.binghamton.edu/zaslav/Tpapers/sggm.pdf)
discusses Harary's balance theorem and cycle obstructions. Our parity-cycle proof
above is included in full and does not assume numerical solver correctness.

Related layered crossing minimization has substantial complexity and
fixed-parameter literature, for example
[Recognizing 2-Layer and Outer k-Planar Graphs](https://arxiv.org/abs/2412.04042).
Results for those graph-drawing variants do not automatically establish hardness
or tractability for this chromosome-block model. No blanket complexity claim or
novelty claim is made here; a publication would need a broader literature review
and independent mathematical review.

## Use

```bash
python certification/structural.py check FIXTURE.json LAYOUT.json CERTIFICATE.json
python certification/structural.py verify FIXTURE.json LAYOUT.json CERTIFICATE.json
python certification/structural.py construct FIXTURE.json NEW_LAYOUT.json CERTIFICATE.json
```

`verify` recomputes the geometric score, bundle counts and cycle witnesses and
accepts only a matching exact certificate. A nonattaining candidate receives
`UNRESOLVED_STRUCTURAL_GAP`. `construct` declines outside the sufficient class.
These operations belong in CPU Slurm jobs on Pegasus for real data. Nothing is
modified in the working optimizer branch, and no new plots are generated.
