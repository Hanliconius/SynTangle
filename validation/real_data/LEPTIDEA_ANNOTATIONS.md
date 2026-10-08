# Leptidea: what annotation preparation means

Target: Fig. 1 of https://doi.org/10.1007/s10577-023-09713-z.
The four male assemblies have downloaded successfully on Pegasus. This does
not establish that their matching gene annotations have been deposited.

## Preferred inputs

The shortest route is the authors' finished MCScanX `.collinearity` and coordinate
files for Fig. 1, together with assembly names/versions and SynVisio settings.
If those are unavailable, obtain the final MAKER gene-coordinate GFF/GFF3 and
matching protein FASTA for each of:

| Assembly | Downloaded accession |
|---|---|
| LsinapisSpaM | GCA_949711715.1 |
| LsinapisSweM | GCA_949711785.1 |
| LjuvernicaM | GCA_949711755.1 |
| LrealiM | GCA_949710795.1 |

The paper also compares Bombyx mori and Melitaea cinxia, with Pieris napi in
additional comparisons. Confirm the exact original comparator versions rather
than substituting current annotations. The methods give an NCBI download date
of 2021-04-25, not a complete accession/version manifest.

## Availability check, 2026-10-08

- [Author repository](https://github.com/EBC-butterfly-genomics-team/Leptidea_chromosome_research2022):
  contains MAKER configuration/merge scripts and BLASTP/MCScanX scripts; final
  annotation and collinearity files were not located in the recursive main tree.
- NCBI Datasets reports for all four versioned accessions returned no
  `annotation_info` and no paired assembly in `assembly_info` during this audit.
  This is evidence about those records, not proof that no annotation exists
  elsewhere.
- The MAKER merge script identifies `.all.maker.gff`, `.all.maker.noseq.gff`
  and merged FASTA outputs. Its input/evidence paths refer to the author's local
  filesystem, so the script alone is not a runnable annotation reconstruction.

No authors have been contacted. A useful request would be the final Fig. 1
collinearity + MCScanX coordinates, or final male-assembly GFF/protein sets,
reference accessions and any scaffold-name conversion table.

## Once those files are obtained

1. Check assembly version and scaffold names against the downloaded FASTAs.
   Archive sequence reports may provide the mapping from deposited accessions
   to original `HiC_scaffold_*` names; require an explicit mapping.
2. Check every protein identifier against the GFF transcript/gene hierarchy;
   record any isoform selection and every discarded record. Check coordinate
   bounds against real sequence lengths. Do not map chromosomes by approximate
   size or silently rename unmatched scaffolds.
3. Generate MCScanX gene coordinates; retain the paper's chromosome selection.
4. Run reciprocal BLASTP, retain the published top-five hits per query, then
   MCScanX with the maximum gene gap of 10. Preserve versions and parameters.
5. Import the resulting observed blocks, recover the paper's species row order
   and chromosome-size display order, then run the comparison and vector PDFs.

Functional annotation with InterProScan/Swiss-Prot is unnecessary for the
crossing objective. We need positions and homologous proteins, not gene names.

## If the original annotations cannot be recovered

Recreating them means repeat masking and multi-round MAKER prediction using
transcript/protein evidence, with SNAP/Augustus training. The paper used MAKER
3.01.04. This is a substantial annotation project and may not recover exactly
the original gene/block set. A new annotation or DNA-alignment route must be
labelled a reanalysis, not reproduction of Fig. 1. Do not start it merely because
the assembly download has finished.
