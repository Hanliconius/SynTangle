# Resume the independently replayed flatworm certificate

Pegasus benchmark 138436391 used a single CPU core:

| Strategy | Exact lower bound | Validated upper bound | Gap | Tree nodes | Search wall time |
|---|---:|---:|---:|---:|---:|
| All cuts, frontier | 55339 | 60294 | 4955 | 101 | 49.571 s |
| Selective cuts, frontier | 58629 | 60294 | 1665 | 101 | 20.123 s |
| Selective cuts, frontier | 59995 | 60294 | 299 | 501 | 76.050 s |

All three certificates passed independent rational replay. The latter remains
unresolved; a 299 gap is not an exact optimality proof. Planarians independently
certified L=U=7901 at the root in the same benchmark.

`submit_flatworm_continuation.sh` snapshots the prior fixture, candidate and
`selective_frontier_501.json` from `proof_combined_BpUBYJ`. It first verifies the
entire saved certificate against its checker fingerprints. `order_resume.py`
then continues its unresolved leaves for up to 900 seconds, with a total tree
cap of 10001 nodes. One compute node, one CPU core, 8 GB, nano, 25 minutes.

The entire old branch tree, selected cut indices, root witness and all exact
node multipliers survive. Closed leaves remain closed. Each new branch partitions
one unresolved domain into both binary cases; both children inherit the parent
witness. Weaker proposed multipliers cannot replace a stronger bound. Floating
infeasibility never closes a branch. The resulting certificate is replayed by
the existing selective checker, and its lower bound must be at least 59995
for the submitted real prior. The original certificate is untouched.

This is actual frontier continuation, unlike the production Leptidea trials,
which start fresh solver searches from a validated best layout. The experiment
remains separate from production. Certificates apply to the stated fixed
block/anchor midpoint-crossing model and permitted whole chromosome moves,
not ribbon-area overlap or every possible rendering convention.
