# Targeted reconstruction repair

Run `bash validation/real_data/submit_simple_repair.sh` from the SynTangle checkout after fetching this branch. It resumes `simple_reconstructed_QI1cVs`, refreshes the comparison renderer, reuses validated completed inputs and pair results, and submits Suzukii and Rhodotorula. It preserves successful Ficus results.

Repeated gene IDs retain distinct placements. Transcript-to-gene resolution requires matching chromosome, strand and coordinate containment; ambiguous parents still fail. Multipart UTR/exon features are not used as gene identity nodes. Original GFF files remain intact.

Rhodotorula now reconstructs all directly observed pairwise alignments for native GENESPACE reference ordering and reference colour projection. Only adjacent-row evidence enters the scored fixture, identical for all methods. Nonadjacent evidence is audited separately. This is reconstructed minimap2 evidence, not the authors' original Syntplot result.

Caméor chromosome selection now accepts the deposited chrNLGN names. Tiffany and Hedin chromosome 1 remains split into L/S in the downloaded source. The targeted retry excludes faba until a documented fixed arm join or equivalent original whole-chromosome representation is available. Independent arm moves or guessed concatenation are forbidden.

Ficus pumila remains source-blocked. Existing successful Ficus comparisons are not rerun by this repair.

The NCBI Dmel source marks mod(mdg4)/CG32491 as trans-splicing with five distinct gene parts. It cannot be represented by one contiguous protein anchor. Its source features remain intact, and the exclusion is recorded explicitly in the input audit. Full-source parsing was checked against GCF_000001215.4, not only a synthetic fixture.

NCBI protein-based inputs now use the supplied protein.faa from the exact assembly's Datasets package rather than translating the whole GFF against chromosome-only FASTA. This avoids requiring excluded unplaced sequences during translation. A new package filename prevents reuse of the previous package lacking proteins. Validation on the actual GCF_000001215.4 protein/GFF package matched 13,985 gene identities from 30,802 supplied protein isoforms before chromosome selection. Existing longest-isoform selection and the explicit trans-splicing exclusion remain in force.
