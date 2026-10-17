# Rhynchospora: a complete small real-data example

This example contains **all 37 regions and 10 chromosomes** from our three-genome GENESPACE reanalysis associated with Figure 1a of [*Sex without crossovers mimics clonal reproduction in Rhynchospora tenuis*](https://doi.org/10.1038/s41586-026-11057-7). No downsampling is used. The two input tables total less than 10 KB; genome downloads are unnecessary.

The rows are *R. tenuis* (2 chromosomes), *R. austrobrasiliensis* (3 selected chromosomes), and *R. breviuscula* (5 chromosomes). Input evidence consists of native GENESPACE `useRegions=TRUE` regions, with fixed directed endpoints and full selected FASTA chromosome lengths.

## Run from the repository directory

Install SynTangle as described in the main README, then:

```bash
syntangle flip-assist \
  --chromosomes examples/rhynchospora/chromosomes.tsv \
  --links examples/rhynchospora/links.tsv \
  --output results/rhynchospora_flips --figure

syntangle run \
  --chromosomes examples/rhynchospora/chromosomes.tsv \
  --links examples/rhynchospora/links.tsv \
  --output results/rhynchospora_full --seconds 20 --figure
```

Both start from the archived GENESPACE layout (13 crossings). Tested release runs reach 9 crossings while retaining all 37 links and 10 chromosomes. Full optimization reports lower bound 9 and gap 0 (`solver_bound_closed`, not an independent exact certificate). The flip helper reaches the same score without supplying a certificate.

To start instead from the selected sequence order (55 crossings), add `--layout examples/rhynchospora/initial_layout.tsv` to either command and choose a fresh output directory. This sequence order is not asserted to be the original paper's exact display order.

From R, without installing GENESPACE:

```r
library(SynTangleR)
x <- syntangle_read_tables("examples/rhynchospora")
flips <- syntangle_flip_assist(x, output = "results/rhynchospora_flips_R")
full <- syntangle_run(x, output = "results/rhynchospora_full_R", seconds = 20)
full$stats
full$pdf
```

## What the figure shows

![Archived Rhynchospora layouts using identical regions](../../docs/images/rhynchospora-example.svg)

The four archived layouts score 55, 13, 9 and 9. The last two are identical here: full search need not outperform flip assistance on every dataset. Panel C flips *R. breviuscula* Chr2_h1 as a whole relative to panel B; it does not reverse individual regions.

This is also a useful illustration of the distinction between crossings and twists. From the GENESPACE baseline, midpoint crossings fall from 13 to 9, while twisted ribbons increase from 11 to 12. That is permitted by the standard midpoint objective. The wide orange ribbon retains its direction in all panels. Optional ribbon refinement optimizes a different objective; see the [user guide](../../docs/user-guide/walkthrough.md).

The documentation figure preserves the archived reference-chromosome colors. The public CLI renderer uses its standard per-link palette, so its colors differ while the evidence and scores agree. Multiple equivalent optimal layouts may exist.

To regenerate the documentation SVG and vector PDF from the archived layouts:

```bash
python examples/rhynchospora/render_example.py
```

This renderer checks each saved crossing score before drawing. It does not rerun the optimization or copy pixels from the published figure.

## Provenance and scope

- Source: `SynTangle_rhynchospora_endpoint_checked.zip`, supplied from the completed Pegasus analysis. The archive hash and table hashes are in `provenance.json`.
- Discovery: three-genome protein-based GENESPACE reanalysis, deposited annotations for *R. tenuis* and *R. austrobrasiliensis*, transferred annotation for *R. breviuscula*. The OrthoFinder compatibility repair preceded these outputs.
- Original source repositories: [paper deposit](https://doi.org/10.17617/3.IXRT5Y) and [companion deposit](https://doi.org/10.17617/3.DGOJRQ). Cite the paper, data repositories and GENESPACE when reusing these derived results. Original source materials retain their respective terms; SynTangle's software license is not a relicensing of those materials.
- This is **not** an exact reproduction of the authors' original twenty-genome discovery. It includes every link in the supplied three-genome region fixture, not every genome or inferred homology from the paper.
- Coordinates, region directions, chromosome lengths and IDs are unchanged. Table conversion encodes a negative block strand by descending endpoints; plotting reflects both endpoints together when a chromosome is flipped.
- `reference_layouts.json` preserves the four saved states and source palette. `reference_optimized_layout.tsv` provides the archived final state. These are reference results, not instructions to rearrange blocks internally.
