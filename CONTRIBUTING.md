# Contributing

Please open an issue describing a reproducible problem or proposed improvement. Include the installed version and the smallest shareable example.

Read [AGENTS.md](AGENTS.md) for the project's development principles: reuse existing code, make the smallest complete change, and retain scientific validation.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test,plotting]'
python -m unittest discover -s tests
R CMD INSTALL r/SynTangleR
SYNTANGLE_TEST_EXECUTABLE="$PWD/.venv/bin/syntangle" Rscript integrations/genespace/test_adapter.R
```

Tests cover legal chromosome moves, scoring, bounds, saved layouts, endpoint preservation, the public CLI and the R hand-off. Build and exercise the installed wheel outside the checkout before changing distribution or entry-point behavior. The CI workflow repeats this on Python 3.11 and 3.12.

## Scientific invariants

Every output must retain the intended chromosome identities, links and internal coordinates. Only whole-chromosome flips and explicitly permitted permutations are legal. Directed endpoints must never be independently sorted during rendering. Do not silently discard links, fill missing native orders, substitute assemblies, or relabel heuristic results as certified optima.

New importers need explicit coordinate conventions and tests for reversed intervals, missing chromosomes, duplicated IDs and malformed records. Preserve provenance and document any filtering or aggregation before optimization.

Keep user documentation in `docs/user-guide`. Historical experiments belong in `docs/development` or validation material; they are not installation instructions. Do not commit virtual environments, caches, private paths, downloaded genomes or newly generated result archives.
