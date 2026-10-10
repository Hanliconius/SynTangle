# Native GENESPACE region comparison

Submit `bash validation/real_data/submit_rhynchospora_native_regions.sh` on Pegasus.
Default source is the existing `local_results/rhynchospora_visual_PTiYgr` run; an
alternate run can be supplied as the first argument. All three jobs use `cpu`.

This addresses the representation mismatch with the author script:

1. Translate deposited Tenuis and Austrobrasiliensis Helixer annotations, and
   the existing transferred Breviuscula annotation, with gffread. Retain one
   longest usable peptide per gene. Preserve the selected 2/3/5 whole chromosomes.
   Audit source hashes and rejected translations. Prefix IDs per species;
   **do not use transferred gene-ID equality to infer orthology**.
2. Run OrthoFinder 2.5.5 and GENESPACE 1.3.1 discovery with MCScanX. GENESPACE's
   within-block orthology and synteny integration remain native.
3. Produce native `useRegions=TRUE` and `useRegions=FALSE` PDFs from the same
   discovery results, explicitly recalculating phased blocks for each. Use the
   author script's five Breviuscula colours, reference order 2/5/1/4/3, and Chr3
   inversion. Save exactly the plotted ribbons, chromosome order and source RDS.
4. Import every adjacent native plotted region once into SynTangle. Require a
   complete chromosome inventory; fail rather than silently drop chromosomes.
   Compare initial, native GENESPACE, native plus our flips, and SynTangle on
   **identical fixed native region evidence**. Checkpoint incumbents and rescore.
   Preserve reversed ribbon endpoints in the interval visualization.
5. Archive both native PDFs, the four-panel comparison, evidence and audits.
   `status.json` distinguishes native discovery completion from solver completion.

The SynTangle PDF uses full selected FASTA chromosome lengths. Native GENESPACE's
own PDF uses its native ribbon-derived lengths and spacing. These display
geometries need not be identical. The scored objective is unweighted region
**midpoint** crossings, not polygon overlap or all original gene-anchor crossings.
Optimality statements apply only to that fixed imported objective.

This is a three-genome **reanalysis**, not a reproduction of the original
20-genome discovery. Breviuscula's annotation is transferred, the precise Fig.1a
chromosome mapping is not verified, and the author's saved region data are not
available. The selected Austrobrasiliensis scaffolds (13,7,1) are supported by
Extended Data Fig.4's pseudohaplotype description. The original author script is:
https://github.com/Raina-M/Rhynchospora_tenuis_project/blob/eae73b27ea4730da6a26477a8edf6dce79f3e011/run_genespace.R

Native homology discovery can take hours; the 300-second cap applies to
SynTangle, after discovery and flip assistance. The separate TE-overlap audit
is not applied in this experiment.

Latest run: `local_results/latest_rhynchospora_native_run.txt`.
Latest archive: `local_results/SynTangle_rhynchospora_native_regions.zip`.
Read `logs/st_rh_native.*.{out,err}` and `logs/st_rh_native_zip.*.{out,err}`.
