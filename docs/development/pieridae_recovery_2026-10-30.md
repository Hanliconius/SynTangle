# Full-evidence Pieridae recovery experiment

The first component audit (job 138662783) found 899 of 912 chromosome
occurrences in one component, including explicitly labelled Z chromosomes.
The twelve frozen components contain eleven singletons and one two-chromosome
strand with optimum 1. They do not establish a detached multispecies Z strand.
The 600-second search improved 14,170,379 to 12,915,397; its lower bound was 1.

## Running the follow-up

Use the existing `validation/real_data/submit_pieridae_components.sh` with the
completed component-test directory as its first argument, an integer search
budget as the second, and `--recover` as the third. This creates a fresh child
run, reuses the installed environment and frozen plotting code, checks the
source receipt/table/layout hashes, and starts from the better complete saved
layout. CPU partition, one core, 16 GB; requested walltime includes one hour
beyond the recovery budget for audits and rendering. Production defaults and
species order are unchanged. No fresh download or general MILP is required.

## Diagnostic outputs

- `labelled_connections.tsv`: all chromosome-pair cells, original labels,
  support counts and fractions relative to each endpoint's partner row.
- `z_incident_connections.tsv`: cells incident to explicitly labelled Z.
  Chromosome 1 is not assumed to be Z.
- `graph_bridges.tsv`: cell edges whose removal disconnects their endpoints.
  These are diagnostic graph separators, not licenses to solve crossing
  optimization independently on either side. Hard orientation couplings and
  ordering interactions require additional boundary states.
- `proposal_chains.json`: chains of reciprocal uniquely strongest partners,
  rejecting any group containing more than one chromosome per species.
  Ties do not choose arbitrary partners. These are not certified homology
  identities; a chain seeded by a Z label can still contain an incorrect match.
- `z_candidate_chain_boundary_links.tsv`: original input rows crossing the
  boundary of the union of Z-seeded candidate chains (when nonempty).
- `recovery.json`: accepted whole-chromosome changes, scores and timings,
  neighborhood termination status and an unconstrained factor lower bound.
- Standard `optimization` receipt, checkpoints and comparison PDFs.

The diagnostic uses the supplied adjacent-row evidence. It does not infer
missing all-pairs homology or assess biological confidence from link counts.
Unlabelled chromosomes outside a Z-seeded candidate chain are not thereby
classified as autosomes. Empty Z reports mean no matching explicit label or
boundary was found, not evidence of biological absence.

## Recovery moves and guarantees

The lossless chromosome-cell tables feed the existing crossing delta cache.
Candidates include coordinated chain placement at either end of its existing
component, coordinated legal chain flips, bidirectional neighboring-row partner
order proposals, legal orientation-basis flips, and whole-chromosome insertion.
Every move retains all original links and coordinates. All accepted moves must
strictly improve the complete objective and pass the separate canonical scorer
and hard-orientation validation. Solved components remain excluded from free
moves, and component order stays consistent. Checkpoints preserve improvements.

Chains remain freely revisable proposals: they are not contracted or frozen.
Neighbor order is not equated with phylogeny or evolutionary events. Species
reordering is deliberately absent because it changes which interfaces are
scored and requires a separate homology reconstruction experiment.

The factor bound is the constant cost plus the sum of independent minima of
all aggregated binary-decision cost tables. Every feasible layout chooses an
entry in each table, so its cost is at least this sum. Inconsistent choices
between tables only weaken this bound. The reported lower bound is the maximum
of that bound, exact presolve contributions and the source's recorded bound
for the identical hashed tables. This is not an independently replayed proof.

A complete no-improvement pass establishes only a local optimum in the listed
neighborhood. A time limit establishes neither a local nor global optimum.
Preparation can itself reach its deadline; that raises an explicit error and
leaves the previously written feasible checkpoint, without a completion claim.
Graph-bridge discovery is not yet an exact separator-state solver.

## Verification

Tiny exhaustive fixtures compare aggregated-cache candidates with canonical
crossing scores and bound them by brute-force optima. Recovery tests verify
monotonicity, hard constraints, move replay, unchanged fixture fingerprints,
and frozen references. Graph tests cover bridge edges and tied partners. The
probe exercises audit, ordinary solve, recovery and the frozen submission worker.
Large Pieridae performance and the real PDF rendering are checked on Pegasus.

## PDF-validation repair

Job 138662941 completed recovery (12,915,397 to 7,835,196; lower bound
334,567; 193 accepted moves) but failed afterward because the nested run did
not inherit the original view's Poppler tools directory. The submission wrapper
now accepts `--validate-figures` as its third argument. It submits a 15-minute,
one-core nano job using the existing frozen plotter's `validate` function;
no optimization or rendering is repeated. Solver files remain untouched.
Poppler is discovered through the resolved tables directory and ancestor views,
or an explicit `SCT_PDF_TOOLS` bin directory. New rendered solve jobs test both
Poppler executables before optimization. Validation-only repairs record their
script hash/commit and resolved executables, and retain original job logs.
