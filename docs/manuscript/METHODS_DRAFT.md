# SynTangle methods draft

Working text for a preprint. Mathematical statements below define the display
problem; empirical claims require the frozen input/result/certificate bundle.

## Fixed evidence and legal layouts

Input comprises an ordered list of species rows, chromosomes of known length,
and fixed homologous block occurrences with genomic intervals. The renderer also
retains the orientation of each occurrence. Homology inference is upstream of
layout optimization. Any filtering or reconstruction of evidence must be described
and frozen before comparing layouts.

For row s, let pi_s be a permutation of its chromosomes and let epsilon_c in
{-1,+1} be the display orientation of chromosome c. A feasible layout comprises
all pi_s and epsilon_c satisfying supplied hard constraints. The legal moves are
whole-chromosome permutations and reversals. Original coordinates, chromosome
lengths and correspondence identities are not optimization variables.

For chromosome c of length L_c, reversal transforms coordinate x to L_c-x.
It acts on every block endpoint on c simultaneously. Directed block endpoints
retain their identities under reflection; sorting endpoint coordinates after a
flip can erase the twist of an inverted correspondence and is therefore forbidden
in directed ribbon rendering.

## Crossing objective

For a fixed layout, let o_c denote the horizontal offset assigned to chromosome c.
For a block midpoint m, its display position is

    X(c,m) = o_c + m                 if epsilon_c = +1,
             o_c + L_c - m          if epsilon_c = -1.

Uniform positive scaling and spacing may be used for display. A homology link h
between adjacent rows has positions X_s(h) and X_(s+1)(h). Two distinct links h,k
cross in the scored objective exactly when

    [X_s(h)-X_s(k)] [X_(s+1)(h)-X_(s+1)(k)] < 0.

The unweighted score sums this indicator over unordered pairs of links in each
adjacent-row comparison. Tied positions do not contribute strict crossings. This
is a midpoint-order objective, not ribbon area, number of ribbon self-twists,
rearrangement distance or an evolutionary reconstruction. A visually twisted wide
ribbon is not by itself a violation or an extra midpoint crossing.

The score applies only to correspondence multiplicities accepted by the canonical
scorer. Ambiguous copies must be resolved under a documented evidence policy or
reported explicitly; they cannot be silently coerced into one-to-one links.

## Optimization and bounds

A feasible state provides an upper bound U on the minimum score C*. A globally
valid relaxation or exhaustive search bound provides L, so L <= C* <= U. An
improved incumbent reduces U without changing the evidence. A better lower bound
raises L. Neighborhood search restricts a local search subproblem and cannot by
itself justify a global lower bound.

The adaptive hybrid path combines legal incumbent improvement with a joint
order/orientation model and bounded neighborhoods. Hard orientation relations
are represented by a GF(2) basis. The default core solver and the experimental
hybrid path must be described separately until their shared release interface is
finalized. Graph/hypergraph structure is useful representation and decomposition
machinery; it does not itself prove optimality.

When the validated global bound equals the independently rescored incumbent, the
solver reports a closed gap for that model. Independent certificate verification
is a separate claim: a checker must validate feasibility, the score, the bound
argument and coverage of the entire legal search space. Floating-point solver
status alone must not be described as a standalone exact arithmetic certificate.
The independent checker is maintained on exact-certification-experiment.

## Comparators and preservation

Baseline and optimized layouts use the same adjacent-row links, lengths and
coordinate conventions. SynTangle flip assistance holds an imported chromosome
order fixed while changing allowed whole-chromosome orientations. It is our
addition, not native GENESPACE behaviour.

The block-proxy benchmark calls GENESPACE's installed internal ordering function;
it is not a complete gene-level GENESPACE rerun. For reference genomes without
qualifying anchors to some chromosomes, that function can omit those chromosomes.
The adapter retains the returned relative order and appends omitted chromosomes
in original input order. Such panels must be labelled supplemented GENESPACE.
All additions are audited and no homology is invented. These results cannot be
reported as unmodified native end-to-end GENESPACE comparisons.

For every rendered panel, checks validate the complete chromosome set, permitted
orientations and the correspondence identity/count. Directed endpoint geometry
is checked against immutable genomic coordinates and whole-chromosome reflection.
The audited PDF is identified by a cryptographic digest. Independent checks of
saved coordinates complement, rather than replace, visual inspection.

## Evidence to freeze before results text

For each case retain accession/version and input checksums, evidence reconstruction
parameters, row order and allowed states, baseline policy, seed and time budget,
actual optimization runtime, feasible layouts, objective scores, global bounds,
endpoint audits and any independent certificate/replay result. Report unsuccessful
or unavailable baselines explicitly. Exact published-figure reproduction and
reconstructed evidence are distinct benchmark scopes.
