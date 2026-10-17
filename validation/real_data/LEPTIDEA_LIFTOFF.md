# Original Leptidea assembly annotation transfer

This pipeline uses the four assemblies from the Chromosome Research paper
https://doi.org/10.1007/s10577-023-09713-z:
LsinapisSpaM (GCA_949711715.1), LsinapisSweM (GCA_949711785.1),
LjuvernicaM (GCA_949711755.1), and LrealiM (GCA_949710795.1).
The annotation donor is the independently annotated L. sinapis assembly
GCF_905404315.1. This is an annotation-transfer reanalysis of the original
assemblies, not a recreation of the authors' original MAKER gene set.

Run `bash validation/real_data/submit_leptidea_liftoff.sh`. An optional first
argument selects the existing `prepared/leptidea` assembly download directory.
The default uses `local_results/latest_published_inputs_run.txt`. Cached target
ZIP files are checked against their original SHA256 manifest and extracted;
target genomes are not downloaded again. The donor FASTA and GFF are obtained
through the NCBI Datasets HTTPS API. Source checksums, commands and installed
tool versions are recorded in the new run directory.

Liftoff 1.6.3 uses coverage >=0.90 and exon/CDS identity >=0.90, excludes partial
mappings and searches for extra copies at identity >=0.95. No chromosome
correspondence list is supplied. All four transfers run concurrently on cpu.
Each worker has a separate annotation database and intermediate directory.

The audit retains donor IDs with one detected target mapping, no detected extra
copies, adequate identity/coverage and nuclear chromosome coordinates. These
are candidate homologies, not independent proof of one-to-one orthology.
Copy discovery is threshold dependent. Unknown IDs, duplicate transfers,
nonchromosomal mappings and low/missing quality values are counted separately.
The initial 90% thresholds are conservative; review the audit before deciding
whether a lower-threshold sensitivity run is needed for divergent assemblies.

Outputs:

- `annotations/ASSEMBLY/lifted.gff` and `unmapped.txt`: full transferred models.
- `homology/annotation_audit.json`: recovery and exclusion counts per assembly.
- `homology/candidate_gene_anchors.tsv`: accepted coordinates, donor IDs and
  mapping quality; coordinates are zero-based, half-open.
- `homology/shared_all_four_gene_ids.txt`: accepted IDs present in all targets.

The TSV includes accepted genes outside the four-way intersection as well.
Collinear block construction and the GENESPACE/SynTangle comparison follow
after reviewing transfer recovery. Neither block boundaries nor the published
layout can be inferred merely from matching donor IDs.
