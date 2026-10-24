# Ponytail principles review — 24 October 2026

## Scope and conclusion

Reviewed SynTangle main at `5c94ea51ee214cdb01bfba2c5efab1609225e865`, Ponytail's README, compact AGENTS instructions and full skill, and the local manuscript/walkthrough/social-figure sources produced in this project. This is a targeted design and validation review, not an exhaustive correctness or performance audit.

We comply well with shared-engine design and preservation requirements, partially with focused verification, and less consistently with reuse and small complete changes in research/figure orchestration. No defensible compliance percentage follows from this review.

Ponytail means minimizing unnecessary work after understanding the real flow. It does not justify removing validation, making mathematical claims weaker or more casual, deleting independently implemented checkers, or replacing readable code with compressed syntax.

## Evidence

| Principle | Assessment | Evidence and consequence |
| --- | --- | --- |
| Reuse existing implementation | Good in public interfaces; weaker in figures | The R package delegates to the shared Python CLI. The local thread and alternative-ending builders duplicate line rendering; the walkthrough adds another version. Reuse rendering semantics on the next edit instead of creating another near-copy. |
| Read callers and fix the root once | Mixed historically | Earlier native-order/schema and endpoint failures required repeated fixes. Current public preservation tests inspect rendered directed endpoints, not only reported crossing counts. Retain that coverage and trace callers before changing adapters or drawing. |
| Small complete change | Mixed | The package has explicit workflows and optional dependency groups, but historical validation and many figure variants remain. File counts alone do not establish bloat: main contains 41 core Python files and 42 top-level test files, and the optimization problem legitimately needs several algorithms. |
| Preserve validation and failure handling | Good public boundaries; one verified research-checker defect | Public workflows reject ambiguous evidence and nonempty output directories. The standalone manuscript checker instead used assertions for required checks; these disappear under Python -O. This defect was reproduced and corrected as described below. |
| Test risky logic | Substantial existing coverage, uneven research coverage | Public tests check illegal evidence, finite budgets, fixed order for flip assistance, preserved links and actual rendered endpoints. Search-pipeline tests compare interrupted bounds/reconstruction against exhaustive small instances. Figure scripts mostly rely on assertions and visual review; these are not a substitute for validating external certificate data. |
| Avoid speculative abstraction/dependencies | Mostly reasonable in the public package | The R interface is thin; plotting and optimization dependencies are optional. Existing JCVI/PAF adapters were explicitly requested, so their presence is not speculative bloat. Their experimental validation limits remain important. |
| Report what is known | Improved, still needs discipline | The public workflow distinguishes solver-bound closure from an independent certificate. Recent figures explicitly distinguish saved comparisons, a fresh replay and an illustrative component composition. A good-looking plot, local convergence or an exhausted budget must not be presented as proof. |

Reviewed implementation paths include `src/syntangle/workflow.py`, `refinement.py`, `public_plot.py`, `visualize.py`, `importers.py`, the R adapter, the endpoint replot script, core CI and the release audit. Reviewed tests include public workflows, release preservation and residual-search continuation.

## Immediate correction: standalone certificate validation

The standalone artifact `source/certify_rhynchospora.py` is distributed with the manuscript and algorithm walkthrough; it is not currently the public package's production solver.

A negative test replaced every saved panel layout with the 55-crossing initial layout while leaving the claimed scores. Ordinary Python rejected the mismatch. Before the fix, Python -O accepted it and printed an exact-certificate status with L=9 and U=55: required assertions had been removed.

Required checks now raise explicit ValueError exceptions in both maintained artifact copies. A focused regression test runs valid and tampered inputs with both normal Python and -O, verifies valid L=U=9 certificates, and requires rejection without a new certificate on tampered inputs. Both copies pass these four cases. The valid original arithmetic result is unchanged. This is a validation fix, not a new mathematical theorem or a change to chromosome rules.

## Targeted follow-up, not a rewrite

1. Consolidate repeated line-drawing and endpoint-transform code when the figure code is next edited. Keep distinct vector/raster backends where useful; share biological coordinate semantics and cover them with endpoint-preservation tests.
2. Give reusable figure builders explicit entry points. Avoid import-time output generation and broad string replacements to maintain variations; retain a small set of named panels rather than another script per presentation.
3. Keep validation snapshots that justify reproducibility, but stop adding new copies of active orchestration. The existing 71,214,796-byte saved-results ZIP is historical material, not grounds for rewriting Git history or deleting evidence.
4. Before removing any method labelled experimental, trace its callers. The public workflow calls the hybrid implementation; filenames are not proof that code is unused.
5. Carry the explicit-validation rule into any publicly shipped certificate interface. Keep its arithmetic checker independent of the optimizer.

Do not rewrite the solver in Rust, replace the optimization architecture, delete the independent checker or add a framework merely to conform to a style preference. Each would need its own concrete benefit and verification plan.

## Persistent guidance and limits

The root AGENTS file adopts these principles while referring to the existing scientific rules. No Ponytail hooks, package runtime dependency or account-wide agent configuration is installed. Future agents must have the repository instructions available to follow them; this change cannot configure unrelated chats globally.

This review did not run the full Python suite, an R/GENESPACE integration, Pegasus workloads, or performance benchmarks. The 17 October release audit records historical checks; it is not evidence that they were rerun today. Production solver and biological model behavior were not changed by this review.
