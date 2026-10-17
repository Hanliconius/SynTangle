# Mandatory directed-endpoint preservation and extended real-data searches

The common comparison renderer now validates legal chromosome permutations and
orientation equations, retains one fixed block evidence signature throughout
rendering, and verifies each individual ribbon against both directed source
endpoints after a whole-chromosome reflection. No endpoints are sorted merely
to make a ribbon look untwisted. Unknown source strand remains explicitly unknown
in the audit. The resulting preservation_audit.json records all fixed link IDs,
chromosome assignments and directed displayed endpoints for each panel, plus the
PDF hash. Completion gates reject missing, partial or stale audits. This is a
layout/rendering legality check, not an optimality certificate or a homology audit.

Both the generic comparison collector and repaired Rhynchospora collector reject
unaudited completed plots. Legacy output needs replotting through the checked
renderer before passing the new gate. The built-in validation submission checks
synthetic endpoint corruption and replots saved Leptidea and Rhynchospora cases;
all three tasks run concurrently, one core each. Full coordinates and evidence
are preserved. Replots do not rerun optimization.

Submit submit_preservation_validation.sh, then submit_leptidea_extended.sh.
The latter automatically depends on successful completion of the former array.
Its three one-core trials run concurrently, with no throttle: adaptive hybrid
900s, global MILP with incumbent hint 900s, adaptive hybrid 1800s. All start from
the same saved best full-evidence state. These are fresh searches, not resumed
proof frontiers. Every incumbent is checked against immutable full evidence.
The collector retains the strongest same-fixture global lower bound and best
legal candidate from both the prior run and all new trials. Conditional or
filtered bounds are never pooled. Result summaries retain solver-reported bounds;
no independent exact certification is implied. If a trial fails, no complete
latest archive is promoted; individual checkpoints remain available.

Logs: st_preserve, st_lep_extend, st_lep_ext_zip. Latest corrected archives are
SynTangle_Leptidea_endpoint_checked.zip and SynTangle_rhynchospora_endpoint_checked.zip.
The extended archive is local_results/SynTangle_Leptidea_extended.zip.
