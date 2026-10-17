# Exact target of the September 22 Bluesky criticism

Post: https://bsky.app/profile/hanliconius.bsky.social/post/3mw4xmy2sgs2q

Paper: Zhang et al. (2026), **Sex without crossovers mimics clonal reproduction in Rhynchospora tenuis**, Nature, https://doi.org/10.1038/s41586-026-11057-7

Target: Figure 1a, not the planarian Figure 4f.

Rows top to bottom: R. tenuis REC haplotype 1 (chromosomes 1,2), R. austrobrasiliensis (1,2,3), R. breviuscula (1,2,3,4,5). Preserve physical scaling and internal inversions. The Bluesky reordered image was a schematic omitting internal inversions; it is not a biological endpoint table.

## Capture the reference

`bash validation/real_data/submit_rhynchospora_capture.sh` submits a cpu job. It preserves the exact published panel as `published_fig1a.pdf`, extracts vector chromosome rectangles in PDF points, archives the Bluesky thread and pinned author script, and inventories the paper and companion deposits. It creates `local_results/SynTangle_rhynchospora_sources.zip`. It does not perform an optimization.

## Sources checked

- Original GENESPACE script: https://github.com/Raina-M/Rhynchospora_tenuis_project/blob/eae73b27ea4730da6a26477a8edf6dce79f3e011/run_genespace.R
- Main deposit: https://doi.org/10.17617/3.IXRT5Y (version 1.0 inspected).
- Companion assembly deposit: https://doi.org/10.17617/3.DGOJRQ.
- Article uses GENESPACE 1.3.1, OrthoFinder 2.5.5, MCScanX 1.0.0 and Helixer 0.3.4 annotations.

No original gene-pair or GENESPACE block table was found in these inventories or the author repository. The published PDF has vector chromosome outlines but raster synteny ribbons; it cannot supply exact biological endpoints. Do not digitize overlapping pixels and call the resulting graph the original data.

## Remaining reproduction work

1. Obtain the exact block output if it exists in another public source. Otherwise rebuild gene-based synteny from deposited assemblies and annotations and label it a reanalysis.
2. Audit the austrobrasiliensis deposited chromosome names and ploidy against the three displayed chromosomes; never silently substitute a different haplotype or all nine chromosomes. The companion deposit has a separately phased hap1 assembly; its identity is not established as the original plot input.
3. Establish the breviuscula assembly/annotation version used in this particular analysis. The companion hap1 FASTA is a candidate, not yet a verified identical input.
4. The published author script runs twenty haplotypes and explicitly changes reference chromosome order and orientation; it is not a turnkey script for the final three-row panel. A three-genome rerun can change inferred orthology relative to the full run and must be labelled accordingly.
5. Render the reconstructed published chromosome order, native GENESPACE, GENESPACE with flip assistance, and SynTangle using the same full evidence. Compare against the archived exact panel before claiming faithful reproduction.

The capture script exposes `--download-assemblies` for a subsequent Pegasus preparation job. It validates deposited file sizes and MD5s and inventories FASTA sequences. This is not enabled by the lightweight capture submission, because assembly identity and annotation completeness remain unresolved. No homology calculation or source capture is run on login nodes.
