# Figure 1a reproduction audit (2026-10-10)

Paper: https://www.nature.com/articles/s41586-026-11057-7
Author script: https://github.com/Raina-M/Rhynchospora_tenuis_project/blob/eae73b27ea4730da6a26477a8edf6dce79f3e011/run_genespace.R

## Confirmed source requirements

The Methods specify Helixer 0.3.4 land-plant annotations, GENESPACE 1.3.1, OrthoFinder 2.5.5, MCScanX 1.0.0, and R 4.2.0. The original discovery involves twenty genomes (two relatives and eighteen tenuis haplotypes). A three-genome discovery is a reanalysis, not an identical orthogroup reconstruction.

The author plotting script uses refGenome=Rbreviuscula, useRegions=TRUE, useOrder=FALSE, braidAlpha=.75, chrExpand=.75, custom palette F9AC60/307BB5/D61F27/ADD8E7/FCF7BF, customRefChrOrder=Chr2_h1,Chr5_h1,Chr1_h1,Chr4_h1,Chr3_h1, and inverts Chr3_h1. These are settings for the all-accession plot; they do not independently establish the final three-row Figure 1a arrangement. The script initializes GENESPACE and plots existing results; it does not include an explicit run_genespace discovery invocation.

## Why current outputs differ

Both new comparisons use Liftoff correspondence instead of orthogroups inferred from independently annotated proteomes. Monotone runs are not GENESPACE reference-phased regions. Two tenuis chromosome colours differ from the breviuscula five-chromosome reference. Keeping all isolated links is not the original region construction rule. Arbitrarily hiding ribbons to match visual density is not an acceptable reconstruction.

## Next implementation requirements

1. Recover original GENESPACE outputs if publicly available; otherwise use deposited Helixer annotations for tenuis and austro, and establish matching breviuscula annotation provenance. Do not substitute transferred donor IDs for inferred orthology.
2. Build proteins and coordinate BEDs with one consistent representative transcript per gene; audit all IDs and coordinates.
3. Execute native GENESPACE discovery (OrthoFinder and MCScanX) and native useRegions=TRUE rendering with the pinned publication settings. Preserve the complete region output returned to the renderer, including exclusions and reference phasing.
4. Verify selected chromosomes, orientations and row order against the published panel. ED Fig. 4 selection alone does not confirm Fig. 1a inputs.
5. Import exactly that retained adjacent-region evidence into SynTangle and rescore every panel on identical evidence. Keep native discovery, ordering-only comparisons and transferred-anchor experiments separately labeled.

The current run-midpoint optimality label pertains only to its reconstructed objective, not the publication evidence. No new native discovery run has yet been submitted.
