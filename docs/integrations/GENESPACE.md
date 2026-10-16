# Refine an existing GENESPACE layout

This is SynTangle assistance for GENESPACE, not a native GENESPACE feature.
Homology discovery is not rerun. Supply a SynTangle fixture JSON containing fixed
block coordinates, chromosome lengths and row order, together with the layout
actually displayed by GENESPACE. Do not treat an arbitrary GENESPACE intermediate
file as a complete display: its schema and orientation conventions must be checked.

The layout TSV has one row per chromosome and these columns:

```text
species chromosome rank orientation
```

Use literal tabs. `species` and `chromosome` must match the fixture IDs exactly.
`rank` is a unique integer within each species, increasing left to right.
`orientation` must be `1` or `-1` relative to the fixture's original coordinates.
Every fixture chromosome must appear exactly once, including chromosomes without
homology. Missing chromosomes, duplicate ranks and illegal orientation assignments
are rejected.

After installing SynTangle:

```bash
syntangle refine-genespace blocks.fixture.json --layout genespace_layout.tsv --mode flips --seconds 60 --output flip_assisted
syntangle refine-genespace blocks.fixture.json --layout genespace_layout.tsv --mode full --seconds 600 --output fully_optimized
```

Flips mode is deterministic greedy improvement, changing only legal whole-chromosome
orientation groups. It preserves chromosome order and reports no optimality proof.
Full mode uses the adaptive hybrid optimizer; it requires the numerical solver
dependencies used by the benchmark environment. Its status is the solver's reported
status, not an independent mathematical certificate.

Both modes export `layout.tsv` and `audit.json`, including the original and final
crossing counts, input hashes, permitted moves and final state. The output is a
layout interchange file: applying it to a GENESPACE plotting object requires an
adapter matching that object's version and coordinate conventions. This command
does not modify a GENESPACE R object or claim to export its native plotting format.
Immutable block coordinates are supplied separately in the fixture; only whole
chromosome order/orientation is returned.

## Unavailable native baselines in benchmark PDFs

If GENESPACE produces no complete layout, the comparison retains all input evidence
and runs SynTangle from the input layout. The two GENESPACE panels explicitly say
unavailable and have null scores/states in `comparison_audit.json`. They are not
counted as wins over GENESPACE. The two rendered evidence panels must pass directed
endpoint-preservation checks. Stale native order files are removed before each run.
