# Real-data starter panel

Source: Lovell et al. (2022), *GENESPACE tracks regions of interest and gene copy
number variation across multiple genomes*, https://doi.org/10.7554/eLife.78526.
The article links https://github.com/jtlovell/GENESPACE_data as its repository of
processed data, intermediate files, scripts, plots and source data.
Pinned archive commit: `77612f8c59fbfd43ef3f4c1719933bf0cbca3261`.

## Selected panel (before looking at optimisation results)

| Table | Selected ordered genomes | Download size | Role |
|---|---|---:|---|
| cotton_split | GbarbadenseA, GdarwiniiA, GtomentosumA | 110 KB | Small subgenome comparison |
| cotton_tetraploid | Gbarbadense, Gdarwinii, Gtomentosum | 171 KB | Intact tetraploid comparison |
| maize | B73, Ki11, Mo18W | 524 KB | Within-species control |
| grasses | rice, Sorghum, maize | 1.31 MB | Between-species plant comparison |
| rho | rice, Sorghum, maize | 2.50 MB | Alternative archived grass block evidence |
| vertebrates | human, mouse, chicken | 6.53 MB | Distant animal comparison |

Sizes are compressed source-table sizes, not predicted search difficulty.
All six tables were inspected for column schema and genome membership. Solver
results and readiness of the selected reciprocal pairs remain pending Pegasus
preparation. These are six cases from ONE paper; grasses/rho also overlap in
biology. They are not six independent studies.

## What requires very little preparation

We can download these block tables directly: no assemblies, protein extraction,
OrthoFinder, MCScanX or all-versus-all alignment reruns are needed. The importer
normalizes the two endpoint spans and emits the existing SynTangle fixture
schema. Every retained directed row becomes one pair-specific observed ribbon.
Reciprocal geometry, orientation and multiplicity must agree exactly before the
reverse-direction copy is excluded. Identical rows within the retained direction
remain distinct links. No transitive orthogroup union is inferred.

The manifest fixes the species chain. Only adjacent pairs are scored; other
pairs/self comparisons are outside the stated benchmark, with source-row counts
and per-edge line/block IDs recorded. Link overlap and repeated chromosome
regions remain represented. Pairwise block orientation is preserved in
provenance rather than treated as a hard chromosome-wide flip constraint.

## What this panel can establish

The first comparison uses the same block evidence for SynTangle and the installed
GENESPACE 1.3.1 chromosome-ordering function, across every reference genome and
two weights. Our flip assistance remains a separate baseline. This is a native
ordering-function comparison on a block projection. It is NOT a reconstructed
GENESPACE gene-level run or a reproduction of the paper's plotted order.

The score counts unweighted crossings of block midpoints, not all constituent
gene pairs or ribbon-border intersections. Initial display order is natural
chromosome-name order. Chromosome extents are maximum observed endpoints,
not full assembly lengths; chromosomes without retained links are outside this
projection. Do not label these initial states "published" or use their extents
as biological chromosome-size estimates. Source strand annotations, gene-scale
homology and assembly lengths need additional import work for full reproduction.

Published biological synteny is already computed, so this panel tests whether
legal display improvements persist on real graphs without mixing that question
with homology inference quality.

## Next datasets to obtain

After this panel, seek independent papers with deposited GENESPACE riparian
objects, gene-coordinate/orthogroup tables or JCVI BED + anchors/simple + seqids
and layout files. Exact published order/orientation is needed for a claim about
excess tangledness of a particular published figure. PDF-only examples are
illustrative leads, not usable benchmark evidence. Include tidy cases alongside
tangled ones; choose the general benchmark on availability/structure before
scoring. Keep a separate illustrative set selected for visible tangle.

Official fallback tutorial: https://github.com/jtlovell/GENESPACE (human/chicken
raw-data example). Processed archive tables above are cheaper starting points.
JCVI format reference: https://github.com/tanghaibao/jcvi/wiki/MCscan-(Python-version)
(requires verified deposited files; not included in this automated panel).

See ../benchmark/FINAL_VALIDATION.md for the Pegasus submission commands and
acceptance gates. No figure generation is performed by either launcher.
