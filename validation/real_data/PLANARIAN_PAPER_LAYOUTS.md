# Planarian paper layouts: panels 4b and 4f

These are separate evidence sets and must not be pooled.

## Fig. 4b: deposited GENESPACE blocks

The first imported fixture contains 585 adjacent block links across 19 chromosomes.
Its original input panel used chromosome-name order, not the published layout.

The labelled Fig. 4b chromosome order is:
- schMedS3h1: 1, 2, 3, 4
- schMedS3h2: 1, 2, 3, 4
- schPol2: 3, 4, 1, 2
- schNov1: 3, 1, 2
- schLug1: 1, 2, 3, 4

The published image places schMedS3h1 at the top, followed by the remaining rows.
Native chromosome coordinate directions are assumed: the figure has no direction arrows
and the author plotting call has no explicit chromosome flip vector. This does not
independently establish every chromosome orientation in the final published image.

Run `bash validation/real_data/submit_planarian_paper_order.sh` after a successful
planarian comparison. It reuses saved GENESPACE, assisted, and SynTangle layouts;
it does not rerun any solver. Output:
`local_results/SynTangle_planarian_paper_order.zip`.

Colours use S. mediterranea h1 chromosome membership. For blocks lacking a direct
reference occurrence, dominant overlap with directly deposited reference blocks
is used for colour only; tied/unassigned mappings are grey. These colour assignments
do not create homology links or change scores. Chromosome lengths remain the imported
maximum observed block endpoints, so this is an order reconstruction, not an exact
reproduction of ribbon geometry or the paper's gene-level display.

## Fig. 4f: deposited BUSCO coordinates

The authors' repository includes all nine full BUSCO tables and chromosome lengths.
No new BUSCO or genome downloads are needed.
Files and plotting rules are pinned to author commit
`91ae879d1d72386d393f38b5e5b4d111a2f1458b`.

The importer follows the chromosome renaming/filtering and keeps exactly one Complete
BUSCO occurrence in every row. An input-only audit found 294 shared BUSCOs,
2,352 adjacent links, and 47 chromosomes, with no empty chromosomes.

The row chain is S. mediterranea h1, S. mediterranea h2, S. polychroa, S. nova,
S. lugubris, Clonorchis sinensis, Schistosoma mansoni, Hymenolepis microstoma,
and Taenia multiceps.

The author script uses numbered chromosome order, equal-width chromosomes,
native directions, and gene-start fractions rounded to three decimal places.
Those display positions define this fixture. Tiny point intervals encode each
position; original base-pair spans and lengths are retained in provenance.
This objective counts crossings of the plotted BUSCO points, not panel-4b block midpoints.

Colours use S. mansoni reference chromosomes with a stable eight-colour palette,
shared across all panels. The GENESPACE baseline uses those BUSCO anchors in the
installed ordering function, not a complete fresh GENESPACE synteny pipeline.

Run `bash validation/real_data/submit_flatworms_fig4f.sh`.
This submits one Pegasus job using the already installed planarian plotting environment,
downloads small pinned input tables, prepares the fixture, runs the baseline/comparison,
and packages `local_results/SynTangle_flatworms_fig4f.zip`.
The solve allowance defaults to 150 seconds; unresolved runs report bounds.

## Sources
- Paper: https://doi.org/10.1038/s41467-024-52380-9
- Fig. 4 image: https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-024-52380-9/MediaObjects/41467_2024_52380_Fig4_HTML.png
- Author Fig. 4b code: https://github.com/Jeremias-Brand/PlanarianGenomeAnalysis/blob/91ae879d1d72386d393f38b5e5b4d111a2f1458b/scripts/genespace_plotting.R
- Author Fig. 4f code: https://github.com/Jeremias-Brand/PlanarianGenomeAnalysis/blob/91ae879d1d72386d393f38b5e5b4d111a2f1458b/scripts/fig4_busco.R

Only importer and syntax checks were performed during implementation.
All solver execution and PDF rendering run on Pegasus.
