# Release candidate audit — 17 October 2026

Audited source: `618e2c359dc71043f5ae6dfe5cc913694c80b36e` on
`release-candidate-0.1`. The accompanying changes add regression checks and CI
coverage; they do not change optimization or plotting behavior.

## Completed checks

- Fresh GitHub clone; isolated Python 3.12 environment; installation with the
  documented optimization and plotting extras; `pip check` passes.
- Built source distribution and wheel. Installed wheel and exercised the CLI
  outside the checkout. Source distribution includes the real example, guide,
  figure, R package, licence and citation; excludes validation results and
  development documentation.
- 162 Python tests pass with development dependencies installed, including
  all four public workflow variants on the real Rhynchospora input.
- 15 additional native-order adapter, reconstruction and preservation tests pass.
- New regression test checks the 74 actual rendered ribbon endpoint tuples
  across input/output panels against the directed source-table coordinates in
  each workflow. It checks unchanged serialized evidence, chromosome/link counts,
  fixed order for flip assistance and the crossing cap for ribbon refinement.
- Installed R package passes schema, original-inversion, table-input and Python
  hand-off checks. Real-data flip assistance and full optimization each give
  13 -> 9 crossings, retaining 37 links and 10 chromosomes, with PDFs.
- `R CMD check`: zero errors, zero warnings, one note for unavailable suggested
  package GENESPACE. R version 4.3.3; Python version 3.12.14.
- README/user-guide/example local links resolve. Licence and citation files are
  present. Public documentation distinguishes numerical bound closure from
  independently verified certificates and marks JCVI/PAF importers experimental.
- Installed CLI PDF visually inspected. Full crossing optimization reports
  13 -> 9 midpoint crossings and 11 -> 12 twisted ribbons: the documentation
  correctly explains why these are separate objectives.

## Remaining release gates

- Run the R adapter against an actual completed GENESPACE 1.3.1 run on Pegasus
  or a colleague's machine. Schema and hand-off tests are not a replacement for
  that end-to-end integration check.
- Review remote CI for the final candidate commit (Python 3.11 and 3.12).
- Independent-user installation smoke test, merge review, version/tag choice,
  release notes and software DOI. No release was published by this audit.

The initial full-test invocation lacked the development-only `cvxpy` dependency.
After installing it, the complete suite passed; no solver assertion failed.
SciPy reports that `threads` is forwarded to HiGHS; this is an expected option
warning, not a failed optimization.
