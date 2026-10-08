# Figure 4f reconstruction

Run `bash validation/real_data/submit_fig4f_reconstruction.sh` on Pegasus.
It reads the completed `flatworms_fig4f_Au9tel` comparison and does not run a
solver, BUSCO, BLAST or homology inference. Optional arguments select the
completed Fig. 4f run and the planarian run containing its plotting venv.

The source is Ivankovic et al., Nature Communications 15, 8215 (2024),
https://doi.org/10.1038/s41467-024-52380-9, Figure 4f. The recipe is pinned to
`scripts/fig4_busco.R` at author commit
`91ae879d1d72386d393f38b5e5b4d111a2f1458b`.

The reconstruction retains the nine rows, native numerical chromosome order,
native orientations, equal-width chromosomes, 0.3 chromosome gaps, rounded
gene-start fractions, chromosome-colour mapping after the R manual colour
scale, thin gene marks and the deposited Bezier control-point equations.
The connecting lines use alpha 0.8 and size 0.5 mm; chromosome bars use 5 mm.
These match the author script's curved-line rendering parameters. The
display-tree topology and orange arrow are reconstructed editorial elements,
not inferred phylogenies or chromosome rearrangements.

The pipeline also downloads the published PDF and archives an original vector
crop of panel 4f. `published_vs_reconstructed.pdf` puts that crop on the left
and the reconstruction on the right. `publication_and_reconstruction.pdf`
contains the original crop, reconstructed native layout, labelled native
layout and saved optimized layout on separate pages. The four `panel_N.pdf`
files render every saved comparison layout with the same author-style recipe.

`reconstruction_audit.json` records source hashes, fixture hash, palette,
saved states, rechecked crossing counts and explicit limits. Native coordinates
are checked against deposited rounded starts. Saved layout validation and
full-evidence rescoring precede rendering. The original gene filter and
provenance remain those of the prepared 294-shared-BUSCO fixture; journal
editing, marker dimensions after figure scaling, and pixel-level equivalence
are not claimed. Inspect the side-by-side PDF before calling the panel an
exact visual reproduction.

Publication figure reused under CC BY 4.0, with attribution and indication of
adaptations in `ATTRIBUTION.txt`. Stable archive:
`local_results/SynTangle_fig4f_reconstruction.zip`.
