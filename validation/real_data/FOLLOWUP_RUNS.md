# Validated-incumbent follow-up

Pegasus batch 138436388 produced Leptidea U=3468 with combined same-evidence
L=359. Preservation batch 138436387 passed for Leptidea and Rhynchospora.
`submit_leptidea_followup.sh` finds the latest completed extended run with
U<=3468 and L>=359, validates its preserved figure, fixture hashes and saved
layouts, then launches two adaptive searches concurrently (1800 and 3600 s).
Both start from the same best saved state. These are fresh searches from a
retained incumbent, not resumed mathematical frontiers. The previous strongest
global bound is retained; neighbourhood or filtered bounds are not substituted.
All evidence and whole-chromosome constraints remain fixed.

The comparison collector retains the best previous/new state and strongest
same-evidence global lower bound, verifies all directed ribbon endpoints, and
writes `local_results/SynTangle_Leptidea_followup.zip`. Incomplete runs do not
promote a complete archive. The prior results are untouched.

The separately implemented flatworm proof continuation is documented on
`exact-certification-experiment`, not integrated into the production solver.
