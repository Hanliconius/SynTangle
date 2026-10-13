# Genome reconstruction comparisons

`bash validation/real_data/submit_simple_reconstruction.sh [INTAKE_RUN]` submits the genome preparation, pair inference, 30/120/600-second optimization comparisons and collector in one operation. Arrays are unthrottled. Each case has independent dependencies. Completed source downloads in the intake run are reused after published-checksum validation.

| Case | Scope | Current reconstruction |
|---|---|---|
| Faba Fig. 2b | Six rows: pea Cameor, pea ZW6, VF8137, Tiffany, Hedin/2, Medicago Mt4 | Supplied annotations and longest protein per gene; reciprocal BLASTP E<1e-5, best subject, MCScanX minimum five anchors |
| Suzukii Fig. 4 | Dmel GCF_000001215.4 and deposited Dsuz; chromosome arms remain fixed units | Supplied CDS/proteins; reciprocal BLASTP E<1e-5, five subjects, MCScanX minimum five anchors; replaces published MMseqs2 |
| Rhodotorula Fig. 3 | Eleven new study strains, in visually verified row order | Minimap2 asm20 primary alignments, MAPQ>=20 and aligned span>=10kb; differs from published Syntplot/BLAST |
| Ficus Fig. 4 | Four independent benjamina comparisons with microcarpa, hispida, hirta, pumila | Same nucleotide reconstruction as yeast; differs from published JCVI protein collinearity |

Ficus pumila CNA0019284 is source-blocked: the inspected CNGB public metadata identifies the assembly but does not expose downloadable genome/annotation assets. Its jobs report SOURCE_BLOCKED; the other three pairs continue independently. No substitute species or assembly is selected. Ficus microcarpa uses GWHABKV00000000 and hispida GWHALOG00000000.1. These versioned reconstructions are not claimed to be the exact files used by the benjamina authors. Original evidence and paper orientations remain unverified.

The earlier intake incorrectly described Ficus comparators and omitted pea ZW6 from faba. The reconstruction manifest corrects those scopes using the actual published panels. Ficus is never converted into an artificial five-row adjacent chain.

Inputs use full selected sequence lengths. Chromosome identities are checked from NCBI reports, GWH FASTA metadata or explicit name patterns and expected counts. Any mismatch fails; no longest-N fallback. Large faba genomes are streamed, rather than collected into chromosome-size strings. Unplaced/excluded sequences are recorded. Gene features retain coordinates and strand; protein identifiers must match the GFF parent graph.

Alignment selection defines the newly reconstructed evidence, before optimization. That evidence remains fixed across all methods. Only whole chromosome permutation and reversal are permitted. Each successful comparison requires the existing directed-endpoint preservation audit. The input panel is labelled reconstructed input, not paper reproduction. A zero solver gap applies only to that fixed objective and evidence.

The collector reports prepared species, completed pairs, audited completed comparisons and source blocks. Its ZIP contains plots, fixtures, manifests and audits, excluding large genome FASTAs. Missing outputs are not counted as complete. Setup installs into a new run-specific environment and does not modify existing environments.

Local validation: five tests cover source scope, refusal to substitute longest sequences, GFF parent identity and fixed coordinates, signed PAF endpoint handling, and sequence-report length rejection. Python compilation and shell syntax checks pass. Downloads, tool installation and full biological inference require the Pegasus run; they have not been executed locally.
