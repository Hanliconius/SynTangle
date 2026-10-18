# Two research directions: species count and element count

Status: elementary proofs under an explicit model, with standard-library
reference algorithms and exhaustive small-case checks. Not peer reviewed or
proof-assistant formalized. Standard dynamic programming, mirror symmetry and
strong-component decomposition are used; their general techniques are not new.
This document extends [the two-row star results](FUSION_FISSION_THEOREMS.md).
It does not claim a complete classification of all two-species instances.

## Model and boundaries

Species rows have a fixed order. Permute and flip whole chromosomes within rows;
keep all supplied link endpoints fixed internally. Only strict midpoint
crossings between adjacent rows are counted; ties are excluded. These theorems
do not optimize ribbon twists, reconstruct evolutionary events, or cover links
scored between nonadjacent rows.

For the species-chain and singleton-separator results, hard orientation
constraints, if present, must be XOR equations within one row or between adjacent
rows. There are no fixed absolute chromosome orientations. Nonlocal hard
equations are explicitly rejected by the reference implementation; they would
require additional states carrying unresolved decisions across a separator.

## Track 1: three, four, then n species

### Theorem 1: an exact recurrence for arbitrary n and bounded row width

For row t let V_t contain every chromosome permutation and flip assignment
satisfying within-row equations. For k_t chromosomes,

    |V_t| <= k_t! 2^(k_t).

For states u in V_t and v in V_(t+1), let c_t(u,v) be their crossing cost if all
adjacent-row equations hold, and infinity otherwise. Define

    D_1(v) = 0
    D_(t+1)(v) = min_(u in V_t) [D_t(u) + c_t(u,v)].

Then the global optimum is min_(v in V_n) D_n(v).

Proof. All constraints and all costs of a prefix ending at row t are represented
by its completed rows and its final row state. Attaching row t+1 adds exactly
one interface cost and the equations on that interface. An optimal prefix ending
in v must have some penultimate state u; its earlier prefix can be replaced by
the optimum D_t(u) without changing any future interface or equation. Conversely,
each finite transition appends a legal row to an attainable prefix. Induction
proves that D_t stores exactly the optimum for every prefix and final state.
Backtracking minimizing predecessors constructs an attaining layout. QED.

This incorporates orientations into the row state. Earlier production layer-DP
code first enumerated complete orientation assignments and then optimized row
permutations. That implementation does not by itself give this new bound as n
grows: the proof experiment now has its own state-local implementation.

There are at most sum_t |V_t||V_(t+1)| transition candidates. For maximum row width
k, this is at most (n-1)(k!2^k)^2, plus state generation and interface scoring.
Direct scoring of an interface with m_t links takes O(m_t^2) comparisons.
For fixed k, the complete algorithm is polynomial in the input size and linear
in species count when per-interface input size is bounded. This is a finite
exact algorithm with an explicit bound, not a branch-and-bound timeout claim.

| Maximum chromosomes per row k | Unconstrained row states k!2^k | Transition pairs per interface |
|---:|---:|---:|
| 1 | 2 | 4 |
| 2 | 8 | 64 |
| 3 | 48 | 2,304 |
| 4 | 384 | 147,456 |
| 5 | 3,840 | 14,745,600 |

These are combinatorial bounds, not measured seconds or a universal practical
cutoff. The reference implementation caps row states and transitions, reports
limits explicitly, and retains prefix tables for exact recurrence review.

### Theorem 2: exactly when adjacent pairwise optima compose

Let p_t = min_(u,v) c_t(u,v), and assume every p_t is finite. Every full layout
has cost at least sum_t p_t. Equality holds if and only if there is a sequence
v_1,...,v_n such that c_t(v_t,v_(t+1)) = p_t for every t.

Proof. Each interface cost exceeds its own p_t by a nonnegative amount. Their
sum is zero exactly when each excess is zero, with the same shared row state
on both of its interfaces. QED.

Thus retain, for each interface, the relation of minimizing state pairs. A path
through these relations is an attainment witness; if no such path exists the
pairwise lower bound cannot be attained. The full recurrence still solves the
problem. Define the compatibility penalty

    delta = global optimum - sum_t p_t >= 0.

This is a precisely defined layout-coupling quantity, not a p-value, biological
rearrangement count or general computational hardness measure.

### Three species: a two-chromosome middle row can obstruct composition

Rows are A (one hub), B (partners a,b), C (one hub). There are two links in each
bundle. Arrange each hub's a and b bundles in separated intervals, a before b.
Give the bundle orientation preferences

    f_A XOR f_a = 0,  f_A XOR f_b = 0,
    f_C XOR f_a = 0,  f_C XOR f_b = 1.

Such evidence is obtained by reversing the internal correspondence in just
the C-b bundle. All four bundle minima are zero. Each adjacent two-species
problem independently has optimum zero, with partner order a,b and suitable
partner flips. Globally XORing the four required equalities gives 0=1, so at
least one bundle must pay one crossing. All-forward flips and partner order a,b
pay exactly that one C-b crossing and no between-bundle cost. The optimum is
therefore exactly 1 and delta=1.

This is an orientation-compatibility obstruction with only two middle
chromosomes. There are no extra hard equations in this example.

### Three species: a purely ordering obstruction

Rows are A (one hub), B (partners a,b,c), C (one hub). Give each partner one
link to each hub, so no bundle has any internal orientation cost. A's positions
are ordered a,b,c, while C's are ordered a,c,b.

A-B has zero crossings precisely for middle order abc or cba, depending on
the A hub flip. B-C has zero crossings precisely for acb or bca. These sets are
disjoint. Thus zero global cost is impossible. Order abc with both hubs forward
costs exactly one crossing on B-C, proving global optimum 1 and delta=1.
Flipping individual middle chromosomes cannot change their relative ranks.

### Four species and arbitrary n

Append a singleton row D and one link C-D to either obstruction. One link
introduces no pair of crossings and imposes no hard equation. The first three
rows still force one crossing, while the same layout extends with no extra cost.
Hence the four-species optimum remains exactly 1 although all independent
interface minima are zero. Repeating this extension gives the same obstruction
for any n >= 3.

This proves that independent pairwise optimality is insufficient once a shared
row imposes incompatible decisions. It does not prove NP-hardness at three or
four species. Theorem 1 still provides an exact algorithm for bounded width.

### Theorem 3: a singleton row is a lossless separator

Suppose row t contains exactly one chromosome and both contiguous subproblems
(rows 1..t and t..n) are feasible. Under the local-equation assumptions above,

    OPT(1..n) = OPT(1..t) + OPT(t..n).

Proof, lower bound. Restrict any full layout to the two subproblems. Interface
costs partition exactly; each restricted cost is at least its own optimum.

Proof, construction. Choose an optimum on each side. Their only shared
chromosome has two possible signs and one possible order. If the signs disagree,
mirror the entire right subproblem: reverse every chromosome order and negate
every chromosome sign. Every strict midpoint comparison reverses on both sides
of every scored interface, preserving crossings and ties. Negating both signs
of any XOR equation preserves its truth. The shared sign now matches. Glue
the two layouts. Locality ensures no omitted constraint crosses the cut. QED.

Repeat the construction at multiple singleton rows. In particular, row sizes

    p, 1, q                 (three species)
    p, 1, 1, q              (four species)
    p, 1, ..., 1, q         (arbitrary species count)

admit exact composition of all adjacent optima. The endpoint star solvers can
handle large p,q using Track 2. This is different from the obstructed row sizes
1,2,1 or 1,3,1 above: it is the structure and location of the shared decisions
that matters, not simply the number of species.

The two-chromosome example makes the singleton guarantee sharp: increasing the
separator from one to two chromosomes can break additivity. It does not mean
every two-chromosome separator causes a penalty.

## Track 2: more partners and more links in low-species problems

The star results remain conditional on compatible preferred bundle flips.
Let w_ij be the exact cost of placing partner i before j on a forward center,
as defined in FUSION_FISSION_THEOREMS.md. Let

    b_ij = min(w_ij,w_ji),  epsilon_ij = |w_ij-w_ji|.

Put an edge i->j of weight epsilon_ij whenever i-before-j is strictly cheaper.
Ties have no preference edge. For any order pi,

    between-bundle cost = sum_{i<j} b_ij
                        + sum_(preference edges backward in pi) epsilon_ij.

Proof. For each pair, its cheaper order pays b_ij; its other order pays exactly
b_ij+epsilon_ij. Summing partitions every between-bundle crossing. QED.

This is precisely a weighted feedback-arc ordering objective on the preference
graph; it may have missing edges because of ties. This connection and the
subset recurrence are standard ordering techniques, not novelty claims.

### Theorem 4: only strongly connected preference groups require cyclic search

Partition the preference graph into strongly connected components (SCCs).
The condensation graph is acyclic. The optimum excess over the pairwise lower
bound equals the sum of the optimal backward-edge penalties inside the SCCs.

Proof, lower bound. Restrict any order to each SCC; its internal penalty is at
least that SCC's optimum, and all remaining penalties are nonnegative.
Proof, attainment. Choose an optimum order within each SCC, and concatenate
SCCs in condensation topological order. Every between-SCC preference edge is
forward, so there is no between-SCC penalty. Ties contribute no penalty. QED.

Apply the exact subset recurrence within each SCC. With d total partners and
SCC sizes s_1,...,s_r, the ordering work after cost preparation is bounded by

    O(d^2 + sum_j s_j^2 2^(s_j)).

For bounded largest cyclic group q this is polynomial as total partner count
grows. Total degree d no longer needs to be small. Acyclic partners have no
cyclic search at all. This improves the earlier implementation that capped the
entire star's partner count. The new cap applies only to the largest SCC.

Example: the 3-partner cycle from the first theorem note plus 97 strictly
ordered partners has 100 partners but largest SCC 3. Its pairwise bound is 12,
and exact optimum remains 13. The recorded subset tables contain 8 entries for
the triple and two entries for each of the 97 singleton SCCs: 202 entries,
rather than a table with 2^100 entries.

More links per chromosome increase cost preparation and integer sizes; they do
not themselves create more chromosome-order states. More chromosomes can create
new ordering states or larger SCCs. These are different axes and should be
measured separately.

## Implementation, evidence and limitations

- `species_chain.py`: full row states with signs and permutations, exact prefix
  recurrence, shared-state compatibility penalty, explicit caps, and optional
  singleton separators. It can dispatch two-row star pieces to the SCC solver.
- `fusion_theorem.py`: preference SCC decomposition, exact subset recurrence per
  SCC, exact geometric rescoring, explicit largest-SCC cap. No numerical optimizer.
- Tests check random three- and four-row fixtures against every legal order and
  flip; both counterexamples; same/adjacent equations; infeasibility; declined
  nonlocal equations; singleton gluing; exact SCC solutions against all orders
  on small weighted matrices; and the 100-partner/3-cycle example.
- The 50-species width-one check uses 196 transition candidates, matching its
  separate bundle lower bound. It does not enumerate 2^50 complete flip states.
- A four-species example with row widths 100,1,1,100 has 202 chromosomes and
  exact score 27. Singleton decomposition gives three proved optima 13+1+13;
  each wide endpoint has only a 3-partner cyclic group. A final mirror alignment
  is necessary and tested. Direct row-state enumeration declines this input.
- General two-row graphs with several high-degree chromosomes are still outside
  the star theorem; small ones are covered by the row-state recurrence.
- General many-species plots with wide rows and large coupled cyclic groups
  remain open for a practically useful structural classification. No universal
  chromosome-number cutoff or new NP-hardness result is claimed.
- Correctness arguments are above; tests validate implementations but do not
  replace the proofs. These computational witnesses are not Lean/Coq proofs.

## Reproduce

Run from the repository root; standard-library Python only:

```bash
python -m unittest discover -s certification -p 'test_species_chain.py' -v
python -m unittest discover -s certification -p 'test_fusion_scaling.py' -v
python certification/species_chain.py certification/examples/species3_order_conflict.json /tmp/species3.json
python certification/species_chain.py certification/examples/species4_order_conflict.json /tmp/species4.json
python certification/fusion_theorem.py certification/examples/fusion_100_partners.json /tmp/fusion100.json --max-cyclic-size 3
python certification/species_chain.py certification/examples/species4_wide_endpoints.json /tmp/species4_wide.json --separators --max-cyclic-size 3
```

## Literature context and next questions

The row-state proof is the standard min-sum dynamic-programming principle on a
chain, applied here with complete chromosome orientation/order states. Its
relation to bounded-width constraint optimization should be positioned alongside
Dechter (1999), [Bucket elimination: a unifying framework for reasoning](https://doi.org/10.1016/S0004-3702(99)00059-4).
The feedback-arc reduction is also standard in crossing minimization; see
Bannach et al. (2024), [PACE Solver Description: UzL Exact Solver for One-Sided Crossing Minimization](https://doi.org/10.4230/LIPIcs.IPEC.2024.28).
Those established general techniques should be cited before any novelty claim.

The next mathematical target is a separator wider than one chromosome where
the number of essential boundary states can still be proved small, or a
multispecies ordered-segment class permitting direct construction. The next
empirical target is to measure row-state counts, largest cyclic group and
compatibility penalty independently of raw link count on controlled examples.
