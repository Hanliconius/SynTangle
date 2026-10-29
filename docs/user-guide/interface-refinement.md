# Experimental interface-guided refinement

Feature branch only; no change to the default optimiser or normative rules.
Built on the conditional companion branch (PR 26). Intended first case: the
Pieris mannii / Colotis evagore interface in the 32-species Pieridae figure.

`syntangle.interface_layout.optimize_interface` accepts a full legal layout and
an adjacent species pair. It reuses the existing pairwise MILP path for 65% of
the requested budget, then propagates its arrangement outward by one-sided
whole-chromosome relocations and orientation moves. Each newly arranged row
anchors the next row. The entire proposal is accepted only for a strict decrease
in the canonical full-family crossing objective. All chromosome copies, links
and native within-chromosome coordinates remain unchanged. Ambiguous homology
is scored as supplied, never converted into a one-to-one chromosome assignment.

This first implementation is sequential outward refinement evaluated as one
coordinated transaction. It does NOT implement rigid chromosome-bundle
contraction or guarantee that the pair optimum can be propagated. Those remain
future experiments. It currently rejects hard orientation constraints explicitly.
The pair solver's original bound and component records are retained; pairwise
closure is not a global certificate. Model construction and canonical scoring
are not preemptible and can exceed the nominal budget. Timeouts are not hardness
proofs. Existing independent running searches and their exclusion ledgers are
never resumed or altered by this companion.

Audit output includes the starting and final layouts, the full proposed layout
(even on rejection), pair solver records, chosen interface, outward row costs,
full proposal score and acceptance decision. The separate Pegasus runner freezes
source hashes and verified checkpoints before beginning and writes the usual
all-pair before/after table. No new proof exclusions are exported.
