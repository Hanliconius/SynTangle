# Exact families with fusion/fission-shaped homology

Status: elementary propositions with proofs under the model below, a reference
construction, and exhaustive small-instance checks. Not peer reviewed or
formalized in a proof assistant. These are applications of standard ordering and
dynamic-programming arguments; no novelty in those general techniques is claimed.
This work belongs to the separate certification experiment, not production.

Continuation: [species-count and element-count extensions](SPECIES_AND_SIZE_THEOREMS.md)
prove a stronger cap based on cyclic preference groups, plus multispecies
recurrences and singleton separators. The implementation's `--max-cyclic-size`
now caps an SCC rather than the total number of partners (`--max-degree` is a
compatibility alias).

## Scope and definitions

Use the exact unweighted strict block-midpoint crossing model in
[STRUCTURAL_THEOREM.md](STRUCTURAL_THEOREM.md). Species rows are fixed; whole
chromosomes may be permuted and flipped. Evidence and internal coordinates stay
fixed. Midpoint ties do not count as crossings. None of the statements optimizes
ribbon twists or infers evolutionary events.

For the main results there are **two species rows**. Collapse all links between
each chromosome pair into one bundle. The resulting bipartite chromosome
incidence graph must be a **star forest**: every nontrivial connected component
has one central chromosome and one or more partner chromosomes. Centers may lie
on either row. Isolated chromosomes are allowed. This includes arbitrarily many
separate fusion/fission-shaped groups, with arbitrarily many links per bundle.
Calling a component fusion-shaped does not establish its evolutionary history.

As before, bundle i has P_i concordant and D_i discordant strict link pairs.
Its minimum internal crossing cost is b_i = min(P_i,D_i). Assume the strict
bundle orientation preferences and any supplied hard XOR equations are jointly
consistent. Without hard equations this is automatic for a star forest.

For a star center c, let A_i be the multiset of native midpoint positions on c
of all links to partner chromosome i. For distinct partners define

    w_ij = number of pairs (x,y) in A_i x A_j with x > y.

Thus w_ij is the between-bundle crossing cost when c is forward and partner i
is placed before j. Partner orientations do not affect this count: every anchor
on the earlier partner chromosome precedes every anchor on the later one.
Multiplicity is retained and midpoint ties are excluded.

## Lemma A: exact reduction of one star to an ordering problem

For a forward center and a partner permutation pi, the minimum crossing cost is

    sum_i b_i + sum_(i before j in pi) w_ij.

For a reversed center, reverse the partner permutation to obtain the same cost.

Proof. Every link pair either belongs to one bundle or connects two different
partners. The first category has cost b_i under preferred relative flips. For
the second, partner chromosome order determines the comparison on one row,
while the center's immutable midpoint order determines the other. This gives
w_ij by definition. Reversing the center reverses every strict comparison;
reversing the partner order also reverses every strict comparison, preserving
crossing status. Ties remain ties. Compatible preferred flips attain all b_i
simultaneously, independently of partner order. These categories partition every
link pair, so no cost is omitted or counted twice. QED.

## Theorem A: a necessary and sufficient condition for pairwise-bound attainment

For every unordered partner pair {i,j}, put a directed edge i -> j exactly when
w_ij < w_ji. Equal costs impose no edge. Call this the preference graph.

Define

    B_star = sum_i b_i + sum_{i<j} min(w_ij,w_ji).

Then an optimal layout attains B_star **if and only if the preference graph is
acyclic**, under the compatibility assumption above.

Proof, lower bound. Each bundle costs at least b_i; each partner pair costs
either w_ij or w_ji, hence at least their minimum. Sum these disjoint categories.

Proof, sufficiency. A topological order exists when the preference graph is
acyclic. It puts every strict cheaper pair in its cheaper order, while ties have
the same cost in either direction. With preferred flips, its cost is B_star.
If the chosen compatible assignment reverses the center, reverse this order.

Proof, necessity. Any total partner order attaining B_star must attain the
minimum of every term: all excess terms are nonnegative. It therefore respects
every strict preference edge. A total order cannot respect a directed cycle.
Consequently a cyclic preference graph makes attainment impossible. Since the
objective is integer valued, a cyclic instance has optimum at least B_star+1.
QED.

The proof does not assert that every acyclic orientation-preference graph in a
general chromosome network suffices. This preference graph is a different graph:
its vertices are the partners of **one star**, and its edges encode ordering
costs, not chromosome flip equations.

## Corollary A1: every compatible two-partner fusion/fission is solved

With at most two partners, the preference graph has at most one edge and cannot
contain a directed cycle. Therefore the optimum has the closed form B_star,
regardless of how many links occur or how interleaved their center coordinates
are. This is strictly broader than requiring separated blocks.

Example. The center's link labels are a,b,a,b; each partner preserves its own
two link positions. Every bundle can achieve zero internal cost. Placing partner
a before b costs one between-bundle crossing; the reverse costs three. Hence
B_star = 1 and the optimum is exactly one, despite bundle-only bound B0 = 0.

## A sharp three-partner counterexample

Take center positions

    A = {2,4,9}, B = {1,6,8}, C = {3,5,7}.

Each partner orders its own three links in increasing center-coordinate order,
so all b_i = 0. The between-partner costs are:

| Order | Cost | Reverse order | Cost |
|---|---:|---|---:|
| A before B | 5 | B before A | 4 |
| B before C | 5 | C before B | 4 |
| C before A | 5 | A before C | 4 |

The preferences are B -> A -> C -> B: a directed cycle. The pairwise bound is
4+4+4 = 12, but no order attains it. Order B,A,C costs 4+5+4 = 13; every total
order violates a strict preference and therefore costs at least 13. Thus the
optimum is exactly 13, proved without trusting a solver or a runtime limit.

This gives an exact boundary for **attaining this bound**: two partners cannot
have a preference cycle; three can. It is not a threshold for computational
intractability. Many three-partner and much larger stars remain acyclic.

## Theorem B: bounded partner count remains exactly solvable with cycles

For a forward star with partner set V, define D[empty]=0 and

    D[S] = min_(v in S) (D[S without v] + sum_(u in S without v) w_uv).

Then the optimum is sum_i b_i + D[V].

Proof. Every permutation of S has a last partner v. Deleting it leaves a
permutation of S without v, and the only deleted crossing terms are w_uv for
the preceding partners u. By induction the best such prefix costs D[S without
v]. Taking the minimum over every possible last partner is exhaustive and
attainable. This proves the recurrence and its optimum by induction on |S|.
Lemma A handles the center orientation. QED.

The direct recurrence uses O(d^2 2^d) arithmetic operations and O(2^d + d^2)
storage, after cost preparation, where d is the partner count. It has a finite,
explicit state bound; there is no branch-and-bound termination assumption.
For fixed maximum d, the class has a polynomial-time exact algorithm as input
size grows. Unbounded d is exponential under this algorithm; that alone is not
a hardness theorem or proof that a better algorithm cannot exist.

A saved recurrence table can be checked with exact integer arithmetic by
checking D[empty] and every displayed recurrence. The reference JSON includes
this table for cyclic stars. This is a simple computational witness under the
proved recurrence, not proof-assistant verification of the implementation.

## Theorem C: combining components is lossless

For two-row star forests satisfying compatible bundle preferences, the global
optimum is the sum of the star optima. Place components in one shared order on
both species rows and use an optimal within-component construction.

Proof. Restricting any legal global layout to one component preserves every
within-component relative anchor order, so its cost is at least that component's
unconstrained optimum. Cross-component crossings are nonnegative. Summation is
therefore a lower bound. Conversely, placing components in the same order on
both rows eliminates every cross-component crossing. Compatible bundle flips
can be chosen globally; whichever orientation they give each center, Lemma A
supplies an optimum by retaining or reversing its partner order. Thus the sum is
attained even when compatible hard equations join different components. QED.

More generally, for any number of rows and arbitrary incidence components, the
same restriction/concatenation argument proves additive decomposition **when
there are no hard orientation equations between different components**. Each
component may contain fusions and fissions and multiple chromosomes per row.
This gives another exact family: bounded total chromosome count per component,
solved by finite enumeration within each component. A coarse enumeration bound
per component of k chromosomes is k! 2^k layouts, times its scoring cost.
The general statement does not turn the two-row star solver into a multirow
solver and does not license arbitrary decompositions of connected components.

## What has and has not been established

- Proven here: a closed-form exact family admitting interleaving fusion links;
  an iff criterion for that bound's attainment; a minimal partner-count boundary
  for possible preference cycles; a cyclic example with an exact positive excess;
  a bounded-degree exact recurrence; lossless component combination.
- Implemented: two-row star-forest recognition, preferred-flip compatibility,
  exact integer costs, acyclic construction, capped subset recurrence, full
  independent geometry rescoring. No production optimizer is imported.
- Not established: a full characterization of tractable chromosome networks,
  novelty of the structural application, hardness at three partners, or a
  multirow fusion theorem. Failure of recognition means only that this theorem
  does not cover the input.
- Next theorem target: multirow fusion/fission components with a shared ordered
  segment structure. Independent optimal ordering of each adjacent-row pair
  cannot simply be assumed to give one compatible global chromosome order.

## Reproduce

From the repository root, standard-library Python only:

```bash
python -m unittest discover -s certification -p 'test_fusion_theorem.py' -v
python certification/fusion_theorem.py certification/examples/fusion_interleaving.json /tmp/fusion_result.json
python certification/fusion_theorem.py certification/examples/fusion_preference_cycle.json /tmp/fusion_cycle_result.json
```

The tests compare constructed scores to enumeration of every legal order and
flip on small examples, including all 90 center words containing two links to
each of three partners. They check both center rows, a reversed center, opposite
direction stars, isolated chromosomes, positive internal costs, the strict
three-partner cycle and explicit rejection outside the hypotheses. These tests
support the implementation; the arguments above are the mathematical proofs.

## Relation to existing theory

The star reduction yields a weighted linear-ordering objective; the subset
recurrence is a standard ordering dynamic program, not a newly invented one.
Classical two-layer crossing minimization already has a substantial complexity
and parameterized-algorithm literature, e.g. Kobayashi & Tamaki (2015),
[doi:10.1007/s00453-014-9872-x](https://doi.org/10.1007/s00453-014-9872-x).
Its ordinary graph vertices have no immutable within-chromosome attachment
coordinates; a reduction is needed before transferring its exact complexity
results to this model. In particular, the star incidence graph alone does not
capture the cyclic preference obstruction above.
