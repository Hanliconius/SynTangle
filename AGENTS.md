# Working on SynTangle

Follow [RULES.md](RULES.md) and [CONTRIBUTING.md](CONTRIBUTING.md). These development instructions do not amend the scientific model.

## Smallest complete change

- Read the affected implementation and trace its callers before editing. Include public Python/R entry points, exports, fixtures, documentation and tests that the change actually reaches.
- Reuse an existing implementation first; then prefer the standard library, a platform feature or an installed dependency. Add a helper, option or dependency only for a concrete need in the requested task.
- Fix a shared root cause once. Avoid parallel scripts for small variations of the same workflow. Keep existing interfaces and conventions unless the task requires changing them.
- Prefer readable code to compressed one-liners. Reduce duplication without weakening validation, error handling, accessibility or scientific auditability.
- Treat complexity as a question to investigate, not a reason to remove working solver paths. Measure a real problem before introducing caches, new algorithms or broad architectural changes.

## Scientific safeguards

- Preserve all chromosome identities, coordinates and supplied evidence under the legal whole-chromosome moves. Share coordinate-transform semantics across renderers; different output backends may remain separate.
- Keep optimizer logic and an intentionally independent certificate checker separate. Independent verification is useful duplication; copied plotting or orchestration logic usually is not.
- Required input and certificate validation must use explicit errors that remain active under Python optimization. Reserve assertions for internal invariants and development self-checks.
- Distinguish feasible layouts, local optima, open bound gaps, numerical solver closure and independently checked exact certificates. A timeout is not a hardness proof.
- Preserve provenance. Label replayed, reconstructed and illustrative material accurately; do not describe them as observed events in an earlier run.

## Verification and scope

- Nontrivial logic needs a focused regression test or meaningful self-check. Test the failure mode and relevant scientific invariants, not a copy of the implementation. Documentation-only edits need no invented test suite.
- When fixing a bug, check every affected caller and report the checks actually run. Protect validation and reconstruction when moving or merging code.
- Keep heavy data preparation, alignment and solver benchmarking on Pegasus. Local work may inspect code and run small unit tests, arithmetic checks and figure rendering.
- Keep one maintained implementation for each supported workflow. Retain reproducibility snapshots deliberately, with their purpose and provenance; do not grow new result directories or archives for every failed attempt.
- Finish with the result, verification and any material untested area. Do not claim checks or capabilities that were not exercised.

Adapted from [Ponytail's development principles](https://github.com/DietrichGebert/ponytail/blob/main/AGENTS.md): choose the smallest change that fully solves the task, while preserving the checks that make it trustworthy.
