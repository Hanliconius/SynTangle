# More proof per search node

This experiment stays on the certification branch. It changes no production
optimizer, evidence, candidate arrangement, or definition of a crossing.

## Two different resources

A **compute node** is a Pegasus machine. A **search node** is a mathematical
subproblem inside one program. The existing proof programs already use one core
per Slurm task. This benchmark packs both real datasets and four variants into
one sequential single-core nano job, avoiding a job array and dependency setup.
It reuses the installed prior proof environment.

## A valid consistency inequality

For binary layout decisions x_u,x_v, write z_uv=x_u x_v. Define

    d_uv = x_u + x_v - 2 z_uv.

On every binary assignment d_uv is XOR. For a four-cycle Q and every odd subset
F of its edges,

    sum_(e in F) d_e - sum_(e in Q minus F) d_e <= |F| - 1.

Proof: XORing edge differences around a cycle gives zero, since each vertex
appears twice. Thus a cycle has an even number of differing edges. Equality to
|F| on the left would require differences on precisely the odd subset F and
agreements everywhere else, which is impossible. Since the left side is integer
and at most |F|, it is at most |F|-1. There are eight inequalities per four-cycle.

Pairwise product relaxations can violate these inequalities. For example, set
all four vertex variables to 1/2; set one edge product to zero and the other
three products to 1/2. All four pairwise McCormick product constraints hold, but
one edge difference is one and the other three are zero, violating the singleton
inequality. Removing this fractional assignment removes no real binary layout.

The interaction graph here joins decisions in adjacent species rows; it is
bipartite. Triangle inequalities would be ineffective, so four-cycles are the
first target. Add only cycles whose four products already exist: **no new
variables**. Rank inspected cycles by objective coupling weights, keep at most
500, and inspect at most 100,000 candidates. This bounded selection is a
performance choice; it has no role in proof validity.

## Preserve bounds and work on the weakest obligation

The previous generator explores depth first and starts children with zero dual
multipliers. The new generator passes the exact parent multiplier vector to both
children. The bound remains valid because a child domain is a subset of the
parent domain. Its box-minimum term cannot decrease when a variable is fixed.
It then selects the weakest remaining lower bound first. This targets the leaves
that constrain the global certificate rather than going deep down one branch.

Every branch still covers both binary values. Neither heuristic arrangements nor
incumbent similarity excludes an ordering. An LP's floating status or reported
objective is never accepted as proof. The checker reconstructs every selected
cycle and replays all multipliers in exact rational arithmetic with residual
corrections. It combines the resulting bound with the independently recomputed
structural bound by taking their maximum, not their sum.

## Controlled Pegasus test

One core, 8 GB, 25-minute nano job. For each saved planarian and flatworm fixture:

- Previous depth-first generator, maximum 101 certificate nodes.
- Weakest-bound-first plus inherited witnesses, maximum 101 nodes, no cuts.
- Four-cycle cuts at the root only, one certificate node.
- Four-cycle cuts plus weakest-bound-first, maximum 101 nodes.

Each generation gets 90 seconds of LP/search allowance; model construction, cut
selection and exact replay are additionally timed. All candidates are identical.
`summary.json` records lower/upper bounds, gap, actual nodes, LP calls, row count,
generation wall time and replay wall time. A stronger root certificate might
replace many search nodes; a larger relaxation might instead be slower. The
experiment must establish which happens. An unresolved result remains unresolved.

Local checks exhaust all 16 binary assignments for the cycle inequalities,
reject a fractional obstruction, test inherited bound monotonicity and reject
certificate tampering. Biological benchmarks run only on Pegasus.

This is a consistency-cut application, not a claim of a new theorem in hypergraph
theory or a general polynomial-time solution to chromosome crossing minimization.
