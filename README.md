# SynTangle

SynTangle reduces crossings in multispecies chromosome-synteny layouts by changing
whole-chromosome display order and orientation. It preserves the supplied homology
links and within-chromosome coordinates.

**Development preview.** The command-line interface and input adapters are being
stabilized for the first release. No stable release has been tagged.

## Install

Python 3.11 or newer is required. From a checkout:

```bash
python -m pip install .
syntangle --help
```

For adaptive hybrid optimization, install the numerical solver extra:

```bash
python -m pip install '.[optimization]'
```

GENESPACE, R, aligners and genome-download tools are needed only for the respective
integration/benchmark workflows, not for the core package or standalone flip mode.
Those workflows have separate environment requirements.

## Try an example

```bash
syntangle analyze-fixture examples/fixtures/layout_noise_only.json --pretty
syntangle optimize-fixture examples/fixtures/layout_noise_only.json --solver local-search --pretty
```

The canonical JSON input records species, chromosomes, full chromosome lengths,
block occurrences and optional hard orientation constraints. See the
[data model](docs/data_model.md), [example fixtures](examples/fixtures/README.md)
and [problem definition](docs/problem_definition.md).

## Refine a GENESPACE layout

Given fixed evidence in the canonical JSON format and a complete exported layout:

```bash
syntangle refine-genespace blocks.fixture.json --layout genespace_layout.tsv --mode flips --seconds 60 --output assisted
syntangle refine-genespace blocks.fixture.json --layout genespace_layout.tsv --mode full --seconds 600 --output optimized
```

Flip mode preserves chromosome order and changes only legal whole-chromosome
orientations. Full mode optimizes both order and orientation. Output includes a
layout TSV, crossing scores and an audit record. The TSV interface is explicit;
a direct GENESPACE R-object adapter is not yet provided. See
[GENESPACE integration](docs/integrations/GENESPACE.md) for the required columns,
coordinate conventions, limitations and treatment of incomplete native layouts.

## What the score means

The current objective counts crossings between corresponding block midpoints on
adjacent species rows. It is not ribbon-overlap area, rearrangement distance or an
inference of the true evolutionary history. Whole-chromosome flips can change
ribbon twists without changing midpoint crossings in the same way.

A feasible layout supplies an upper bound. A valid global lower bound supplies a
limit on possible improvement. Matching solver bounds are reported separately
from independent certificate verification. Local/neighborhood optimality alone
does not establish global optimality. Claims apply only to the fixed input,
objective and allowed display states.

No optimizer may split a chromosome, rearrange internal blocks or delete homology
to improve the score. [RULES.md](RULES.md) defines these invariants. Directed
endpoint-preservation tests protect the plotting transformations.

## Reproducibility and development

- [Validation and benchmark workflows](validation/benchmark/README.md)
- [GENESPACE integration and baseline policy](docs/integrations/GENESPACE.md)
- [Release preparation and evidence status](docs/RELEASE_PREPARATION.md)
- [Historical development overview](docs/development/README_before_release_cleanup.md)
- [Benchmark history](docs/benchmark_history.md)

```bash
python -m pip install '.[test]'
python -m unittest discover -s tests
```

Large genome downloads and biological benchmark jobs are designed for an HPC
scheduler. Reconstructed homology comparisons are labelled separately from exact
reproductions of published evidence. Missing baseline results are not scored as
wins. Experimental independent certification remains on the separate
`exact-certification-experiment` branch.

## License and citation

SynTangle code is distributed under the [MIT license](LICENSE). Third-party
software and source datasets retain their own terms. A versioned archival
DOI/citation is a release prerequisite; no DOI is claimed here.
