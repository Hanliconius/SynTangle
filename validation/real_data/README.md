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
Reciprocal geometry and multiplicity must agree exactly before the reverse-
direction copy is excluded. Both directional orientation annotations are retained;
disagreements are flagged explicitly and never become hard flip constraints. Identical rows within the retained direction
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

## Successful final Pegasus validation and the first comparison PDFs

User-reported run 138413416 passed all three groups: 144 discovery tests,
2 exhaustive tests, and 35 focused tests (groups overlap). Run 138413417 imported
all six panels, including rho with one explicitly preserved orientation conflict.
This establishes preparation and validation, not real-panel solver results.

Submit comparison and vector plotting on Pegasus only:

```bash
bash validation/real_data/submit_comparison.sh
```

The launcher defaults to `local_results/latest_final_run.txt`; an explicit final
run directory can be supplied as its first argument. It reuses prepared inputs,
snapshots code/data, installs plotting dependencies in a separate venv, then
submits all six cases concurrently. Default solve allowance is 150 seconds per
case, configurable with `SCT_REAL_SECONDS`. GENESPACE 1.3.1 is read from
`lep_busco_painter_clean`, configurable with `SCT_REAL_GS_ENV`.

Each case writes `comparison.pdf`, saved states, bounds and audit JSON. Four
panels show natural-name input, the installed GENESPACE ordering function on
block proxy anchors, our fixed-order flip assistance, and SynTangle. Input and
baseline PDFs are checkpointed before solving. An interrupted job may therefore
leave a partial comparison PDF; only a case with `COMPLETE` has finished all four
panels. The launcher does not recreate the paper's plotted order or gene-level
GENESPACE pipeline. All source links and overlapping block spans are retained.

The block projection contains pair-specific links rather than multispecies
orthogroups. End-species references may therefore return incomplete orders.
The real-data launcher skips such variants wholesale and records missing
chromosomes/errors in `native_variant_audit.tsv`; it never fills an omitted
chromosome with an invented position or drops it from the fixture. At least one
complete native variant is required. Labels and audit identify this restricted
comparison. The synthetic benchmark helper remains strict by default.

An `afterany` collector creates `SynTangle_real_comparisons.zip` with individual
PDFs, prepared inputs, saved layouts/results/audits and an explicit completion
manifest. Its combined PDF includes completed cases only. The ZIP is also copied
to `local_results/SynTangle_real_comparisons.zip` for a single stable download
path. Download from a local terminal using the same host/alias as normal SSH:

```bash
scp 'jjhanly@YOUR_PEGASUS_HOST:/absolute/checkout/path/local_results/SynTangle_real_comparisons.zip' .
```

See [PAPER_DATA_ACCESS.md](PAPER_DATA_ACCESS.md) for access and preparation status
of independent published figure candidates, including the required 2023
Leptidea example.

## First Pegasus preparation (final_NFyQeq)

Five panels imported: cotton_split 71 links / 39 linked chromosomes; intact
cotton 663 / 80; maize 36 / 30; grasses 453 / 32; vertebrates 871 / 80.
The rho importer initially rejected Sorghum/maize reciprocal annotations.
Inspection of the pinned source found 618 links in each direction with identical
geometric multisets, but one link has + versus - orientation annotations:
Sorghum Chr08:56666140–56906563 / maize 4:1332953–3505204.
The updated importer preserves both annotations and flags this disagreement.
Its block-midpoint objective does not use those annotations as hard parity.
Geometric or multiplicity discrepancies still cause import failure. This is
a source annotation discrepancy, not an optimisation result.

## Published-paper scope

See [PAPER_DATA_ACCESS.md](PAPER_DATA_ACCESS.md) for the active shortlist.
Use `bash validation/real_data/submit_published_inputs.sh` on Pegasus to prepare
the deposited planarian block projection and download the four Leptidea male
assemblies concurrently. Leptidea still needs gene annotation/collinearity
reconstruction; the planarian import is not yet an exact published-display
replica. Other candidates are deferred until finished plotting inputs suffice.

## Planarian comparison after preparation

```bash
bash validation/real_data/submit_planarian_comparison.sh
cat logs/st_planarian_fig*.{out,err}
```

The script reads `local_results/latest_published_inputs_run.txt`, submits setup,
comparison and collector jobs, and reuses the prepared data. An optional first
argument chooses a specific published-input preparation run. SynTangle gets
150 seconds by default (`SCT_REAL_SECONDS`); Slurm allows 20 minutes for the
complete comparison, baseline trials and PDFs. Four vector panels compare input
chromosome-name order, GENESPACE ordering, GENESPACE plus our flips and SynTangle.
All 585 scored links are retained. Deposited nonadjacent blocks supplement
GENESPACE's reference ordering only; their use is logged in
`native_input/reference_anchor_audit.json`. This is the installed GENESPACE
1.3.1 ordering function on block proxy anchors, not the original 1.0.8
gene-level plotting run. Original published chromosome order/flips remain
unrecovered. The collector bundles results, bounds, provenance and PDFs at
`local_results/SynTangle_planarian_comparisons.zip`; inspect its completion status
before interpreting partial PDFs.

For the remaining Leptidea work, see [LEPTIDEA_ANNOTATIONS.md](LEPTIDEA_ANNOTATIONS.md).
