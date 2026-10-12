# Codfish and Carex reanalyses

`submit_next_round.sh` submits these new cases, the Leptidea follow-up and
the independently cloned flatworm certification continuation. No search
executes on the login node.

User-selected cases: Genome Biology DOI 10.1186/s13059-026-03975-6, Fig.2A,
and Molecular Ecology DOI 10.1111/mec.17086, Fig.3. PNAS DOI
10.1073/pnas.2410939122 is removed from the active set.

`submit_codfish_carex.sh` submits independent preparation arrays for the six
codfish and four Carex-paper species, then 15 and six reciprocal protein-search
jobs, respectively. Arrays have no throttle. Each dataset waits only for its own
preparation. This pipeline translates supplied gene annotations, chooses the
longest usable translated protein per gene, runs reciprocal BLASTP and MCScanX
(maximum gene gap 10), and compares fixed signed blocks using the installed
GENESPACE ordering function, our flip assistance and adaptive SynTangle.
The resulting GS comparison is ordering on block proxy anchors, not a full
native gene-level GENESPACE discovery workflow.

Sources:
- https://zenodo.org/records/17475568: codfish FASTA and GFF, pinned deposited MD5s.
  NEAC uses GCF_902167405.1 FASTA and the deposited NEAC annotation.
- https://zenodo.org/records/6990077: Carex cristatella and scoparia final FASTA
  and GFF, pinned deposited MD5s. These have 35 and 33 assembled chromosomes;
  33 is the assembly count, not the observed 31-chromosome cytological count.
- Carex littledalei: resolve the assembly linked to deposited CM022079, retain
  the exact published CM022079-CM022107 chromosome accession set, and require
  an annotation in the same versioned NCBI package. No species substitution.
- Rhynchospora breviuscula: GCA_027562975.1; require all five assembled nuclear
  chromosomes and an annotation from the same package.

A parallel inventory downloads metadata from the author's three analysis
archives (Zenodo 8138413/8138415/8138417), the Carex annotation deposit and
littledalei Figshare 12197544 v1. Metadata only, not indiscriminate archive
extraction. It records actual available filenames for annotation recovery if
NCBI lacks the requested GFF. A missing GFF or ambiguous sequence mapping stops
that species' preparation; its downstream jobs cannot produce a complete plot.
Do not interpret a blocked Carex pipeline as an optimization failure.

The codfish Fig.2A row chain was inspected: Atlantic cod NEAC, haddock, polar
cod, Arctic cod, burbot, hake. Chromosome sizes determine ordering outside
Atlantic cod, matching the caption. Atlantic cod's displayed permutation and
whole chromosome orientations are not recovered yet. The Carex four-species
chain is a study-based chain; exact figure permutations, row direction and
orientations remain to be checked. Both inputs are clearly labelled reanalyses,
not replicas of the original drawings. New blocks cannot establish the original
paper's crossing count. Retain published-layout comparison as a separate goal.

Fixed evidence: full selected chromosome lengths; all adjacent blocks and signed
endpoints; all selected chromosomes including unlinked ones. Every exclusion
of unplaced sequences is recorded before inference. No smallest-hit filtering
is performed here. If native ordering omits chromosomes, an incomplete variant
is rejected, never silently completed or used to claim a win. The shared
preservation audit must pass before a comparison is marked complete.

Read small summaries in `logs/st_codfish_prep.*`, `logs/st_carex_prep.*`,
`logs/st_carex_sources.*`, then corresponding `*_fig.*` logs. Results are isolated
under `local_results/codfish_carex_*`. Successful bundles are
`local_results/SynTangle_codfish_comparisons.zip` and
`local_results/SynTangle_carex_comparisons.zip`.

Neither selected study is assumed to be computationally simple. The conserved
Carex scoparia/cristatella pair is a potential easier control, but the main
four-species figure retains the other paper species.
