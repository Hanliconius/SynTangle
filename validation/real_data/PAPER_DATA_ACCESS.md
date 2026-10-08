# Published figure data: access audit, 2026-10-08

Preparation difficulty is an estimate, not solver difficulty. A downloadable
genome does not establish availability of the exact homology links drawn in a
paper. For an exact before/after demonstration we need those links, endpoint
coordinates, chromosome lengths, species row order, and displayed chromosome
order/orientations. Recomputed links are a separate reanalysis.

## Priority candidates

| Study / figure | Verified access | Remaining work | Estimated preparation |
|---|---|---|---|
| [Leptidea, Chromosome Research 2023](https://doi.org/10.1007/s10577-023-09713-z), Fig. 1 | Paper explicitly orders chromosomes by size. [Author repository](https://github.com/EBC-butterfly-genomics-team/Leptidea_chromosome_research2022) has BLASTP filtering, MCScanX and plotting scripts. ENA PRJEB58697 returns all eight assembled genomes. | Finished collinearity files and the gene annotations required for Fig. 1 were not located in the repository tree. Obtain `.collinearity`, MCScanX coordinate inputs, chromosome lengths and SynVisio settings from authors; alternatively verify annotations and rerun the published pipeline (maximum gene gap 10). Linkage-map files are not a substitute for this figure's links. | Moderate if processed files supplied; substantial if rebuilding. Required illustrative case. |
| [Planarians, Nature Communications 2024](https://doi.org/10.1038/s41467-024-52380-9), Fig. 4b | Paper explicitly supplies inferred GENESPACE blocks in Supplementary Data 13 and 14, within the Supplementary Data 1–15 workbook. | Inspect table schemas, import the selected pairwise blocks, retrieve chromosome lengths, and reproduce displayed order/orientation. Workbook not yet converted or validated. | Low to moderate; strongest independent ready-block lead. |
| [Rhynchospora tenuis, Nature 2026](https://doi.org/10.1038/s41586-026-11057-7), Fig. 1a | [Repository](https://github.com/Raina-M/Rhynchospora_tenuis_project) identifies `run_genespace.R` as the Fig. 1a script. It specifies GENESPACE settings, a custom reference chromosome order and one manual inversion. | Repository contains scripts but no finished block tables/GENESPACE object in its tree. Request the plotting object, block table, BED and chromosome lengths; confirm exact final panel order, since the script also plots all accessions. | Low with saved objects; moderate/substantial rebuild. Ideal 2/3/5-chromosome illustration. |
| [Codfish, Genome Biology 2026](https://doi.org/10.1186/s13059-026-03975-6), Fig. 2 | [Zenodo 17475568](https://doi.org/10.5281/zenodo.17475568) file list verified: assemblies, GFFs and protein FASTAs, including gadMor3/NEAC annotation. | No finished MCScanX collinearity files in archive listing. Run published protein homology/MCScanX workflow or request finished inputs; reproduce plotted ordering and filtering. | Moderate; strong non-Lep complex candidate. |
| [Carex, Molecular Ecology](https://doi.org/10.1111/mec.17086), Fig. 3 | Zenodo [8138413](https://doi.org/10.5281/zenodo.8138413), [8138415](https://doi.org/10.5281/zenodo.8138415), [8138417](https://doi.org/10.5281/zenodo.8138417) each archive a species analysis repository. Public repository trees contain R notebooks and linkage/feature data. | No finished GENESPACE synteny object/block file located in the inspected trees. Do not confuse linkage-map plotting inputs with Fig. 3 inputs. Request original GENESPACE outputs and exact assembly versions. | Uncertain, provisionally moderate/substantial. |
| [Biscutella, Nature Communications 2026](https://doi.org/10.1038/s41467-026-73793-8), Fig. 1d | [Figshare 28451828](https://doi.org/10.6084/m9.figshare.28451828) verified: eight assemblies and GFFs. Paper puts ancestral syntenic fragment definitions in Supplementary Data 5; uses NGenomeSyn. | Inspect fragment table for coordinates and copy identities. Ancestral block labels alone must not be merged into ambiguous pairwise links. Prefer original NGenomeSyn link/layout files. | Low/moderate if fragment table is sufficient; otherwise moderate rebuild. |
| [Melinaea/Mechanitis, PNAS 2025](https://doi.org/10.1073/pnas.2410939122), Fig. 4A | [Author repository](https://github.com/rapidspeciation/mechanitis_melinaea) contains synteny plotting scripts, including PAF refinement; assembly accessions reported in paper. | Inspected repository tree did not expose finished PAFs. Obtain filtered plotting PAFs or align exact assembly versions and apply published filtering. | Moderate; no need to repeat assembly or raw-read analyses. |
| [Codonopsis, Frontiers Plant Science 2024](https://doi.org/10.3389/fpls.2024.1469375), Fig. 2B | [Figshare 26799025](https://doi.org/10.6084/m9.figshare.26799025) verified: assembly, GFF3, CDS and protein FASTA. | Finished comparative links absent from this archive listing. Identify exact other-species annotations and reproduce the published synteny workflow or obtain link tables. | Moderate. |
| [Cordyceps, Scientific Data 2026](https://doi.org/10.1038/s41597-026-07231-1), Fig. 8 | Paper deposits focal assembly/annotation at CNCB GWHGQKB00000000.1 and assembly at NCBI JBRIOS000000000.1. | Finished multispecies link files not verified; locate comparator versions, synteny parameters and links. | Provisionally moderate; metadata availability verified, download/import not yet tested. |
| [Cryptococcus/Kwoniella, PLOS Biology 2024](https://doi.org/10.1371/journal.pbio.3002682), Fig. 1C | [Zenodo 11199354](https://doi.org/10.5281/zenodo.11199354) archive and corresponding GitHub tree verified. Includes orthogroup tables and several figure scripts. | README lists Fig. 1 phylogeny, but not the Fig. 1C synteny rendering inputs. Orthogroups alone lack ordered endpoints. Obtain block/layout files or pair with exact annotation coordinates and rebuild. Existing manual layout improvements make it a useful tidy/control candidate. | Moderate/substantial for exact replication. |

## Earlier exploratory leads: not yet qualified as ready data

These remain figure leads, not accepted reproducible inputs. No claim that their
final plotted links have been located:

- [Cucumis hystrix](https://doi.org/10.1038/s41438-021-00475-5), Fig. 3a:
  paper describes MUMmer/MCScanX; inspect Supplementary Tables S7–S9 for actual
  coordinate coverage before calling them ready. Provisionally low/moderate.
- [Muntjac deer](https://doi.org/10.1038/s41467-021-27091-0), Fig. 1d:
  author code at https://github.com/YinYuan-001/muntjac_code and
  https://doi.org/10.5281/zenodo.5533097; exact synteny inputs still unverified.
- [Sugarcane](https://doi.org/10.1038/s41586-024-07231-4), Fig. 2c:
  defer until copy-specific links and haplotype handling are verified. Complex
  polyploid evidence is a distinct import problem; not a minimal-prep first case.
- [Leptidea 3D chromatin](https://doi.org/10.1038/s41467-024-50529-0), Fig. 1e:
  compare assembly versions and ascertain whether the 2023 processed inputs can
  be reused before creating a second apparently independent dataset.
- [Parotis](https://doi.org/10.1038/s41597-025-05053-1),
  [beet webworm](https://doi.org/10.1038/s41597-025-04371-8), and
  [Clanis](https://doi.org/10.1038/s41597-024-03853-5): assembly/data descriptors,
  but exact plotted links not yet located. Clanis annotation deposit:
  https://doi.org/10.6084/m9.figshare.25151900.v1. Provisionally moderate rebuilds.
- [Cephalopods](https://doi.org/10.1038/s41467-022-29748-w): retain as a possible
  graph-complexity dataset; earlier inspected panel was a dotplot, not the
  riparian before/after target. Import readiness unverified.

## Immediate order of work

1. Plot/solve the six already prepared GENESPACE paper block projections.
2. Import planarian supplied blocks and inspect Biscutella fragment coordinates.
3. Retain Leptidea as a required target; secure its exact Fig. 1 inputs.
4. Secure Rhynchospora Fig. 1a objects for the small illustrative demonstration.
5. Rebuild codfish links if originals are unavailable.

No authors were contacted during this audit. No real-data solver or figure
generation was run locally.
