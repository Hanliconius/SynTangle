# Rhynchospora Fig. 1a: input recovery

Run `bash validation/real_data/submit_rhynchospora_inputs.sh` in an activated
Python environment on Pegasus. All five download tasks run on cpu without an
array throttle. The collector runs after all tasks, including failures, and
prints a short inventory. A failed download makes the collector fail too.

This is input recovery, **not an optimized comparison**. It downloads the REC
hap1 FASTA/GFF, austrobrasiliensis FASTA/GFF, and companion breviuscula hap1
FASTA from Edmond. Sizes and MD5 checksums are pinned from the inspected
deposit inventory. Existing files are checked too; invalid caches fail rather
than being accepted. An optional run-directory argument reuses those downloads.
Downloads allow up to three one-hour attempts and resume retained partial
files after connection/time-out failures. Only verified files are promoted
to the final filename. A server that refuses range requests triggers a restart
of that partial; failed transfers otherwise remain available for a later retry.
The download jobs have a four-hour wall-time allocation, not extra memory.
Completed files are verified and reused when the same run directory is supplied.
No large source archive is produced. Detailed names/lengths and GFF coordinate
checks are saved in `input_readiness.json`.

Three uncertainties must remain explicit:

- The published panel has 2, 3 and 5 bars. Deposited sequence counts are not
  sufficient to establish their identity or order. No scaffolds are joined,
  excluded, or renamed by this stage.
- The author's supplementary Fig. 14 script names multiple austrobrasiliensis
  scaffolds within chromosome groups. Their relationship to Fig. 1a needs
  validation before treating them as whole-chromosome movable units.
- The inspected deposits did not supply a matching breviuscula GFF or the
  original GENESPACE block table. A separately deposited older reference GFF
  must not be attached directly to the new assembly. Annotation transfer or
  an explicitly described genome-alignment analysis is needed if the original
  inputs cannot be recovered. Recomputed evidence is not an exact replica of
  the published ribbons.

Source identities: Nature DOI 10.1038/s41586-026-11057-7; Edmond deposits
10.17617/3.IXRT5Y and 10.17617/3.DGOJRQ. Author scripts inspected at
`Raina-M/Rhynchospora_tenuis_project` commit
`eae73b27ea4730da6a26477a8edf6dce79f3e011`.

The production optimizer and independent certification branch are unchanged.

## Annotation and mapping follow-up

`bash validation/real_data/submit_rhynchospora_followup.sh` uses the verified
inputs from `rhynchospora_inputs_UkAKU9` and geometry from
`rhynchospora_sources_ACs9FQ`; optional arguments override both directories.
Setup, Liftoff, and an annotation audit run on Pegasus CPU nodes. Detailed
environment and Liftoff logs stay in the result directory. This code has syntax
checks only here; biological execution and validation happen on Pegasus.

The exact breviuscula assembly receives a **new** annotation transferred from
the deposited tenuis GFF. Liftoff uses 0.5 coverage/identity to record a broad
candidate set, without its extra-copy search; a separately retained strict set
requires unique recognized donor gene IDs, at least 0.9 coverage and identity,
no partial-mapping flag, and valid target coordinates. This subset is not an
orthology proof and does not recreate the publication's gene predictions.
All transferred GFF records remain available alongside the explicit exclusions.

The mapping diagnostic ranks ordered three-sequence selections from the nine
austrobrasiliensis sequences by similarity to the published bar-length ratios.
It writes ten candidate mappings and the entire original sequence inventory.
Length matching is suggestive only; none of these candidates is automatically
selected, renamed, joined, or passed to an optimizer. A different plotted scale
or a composite chromosome representation can invalidate this diagnostic.

No optimized PDF is claimed by this stage. The original GENESPACE block table
and confirmed chromosome mapping are still required for exact figure replication;
an analysis of recomputed evidence must be identified as such.

## Reconstructed three-species comparison

The paper's Extended Data Fig. 4 caption explicitly states that scaffolds 1, 7,
and 13 were used as an austrobrasiliensis pseudo-haplotype because its six
haplotypes were not fully resolved. This establishes a cited selection for that
analysis, not an exact identification of Fig. 1a inputs. The corresponding
`a1` sequences are annotated. Their 13,7,1 starting order is suggested by the
short-long-long bar widths; orientations remain native and unverified.

`bash validation/real_data/submit_rhynchospora_comparison.sh` creates a new
snapshotted run with four dependent nano jobs, all within the 30-minute partition
limit. It reuses `rhynchospora_followup_jY8c5Y` by default (optional argument
changes that path), including the completed breviuscula transfer and Liftoff
installation. It extracts the three selected austro sequences without modifying
the deposit or indexing it in place, and transfers the exact same tenuis donor
annotation to this selected representation. A timed-out transfer must be rerun
with an appropriate wall time before any comparison can be accepted.

Both transfers use uniquely recognized donor IDs, valid target intervals,
coverage and identity >=0.9, and no partial-mapping flag. All strict Tenuis-Austro
anchors and all strict shared-donor Austro-Brevi anchors are retained. Brevi-only
anchors have no adjacent Austro endpoint in the fixed row chain and are counted
explicitly as nonadjacent exclusions. Shared IDs mean transfer correspondence,
not independently validated reciprocal orthology. The original GFFs and explicit
quality exclusions are retained; these links do not recreate published GENESPACE
orthogroups or block boundaries.

The full 2/3/5 selected chromosome representations and exact FASTA lengths are
preserved, including chromosomes lacking scored anchors. The objective is strict,
unweighted gene-midpoint crossings of adjacent species. GENESPACE ordering
functions, the additional flip heuristic, and SynTangle use identical retained
anchors and are rescored on this objective. The starting layout is labeled a
candidate paper order, never the published drawing. Colours follow the two
Tenuis donor chromosomes and remain stable across all panels. This is an anchor
ordering comparison, not a rerun of the original GENESPACE discovery pipeline.

The comparison uses a 300-second hybrid solve and saves honest bounds. Completed
PDFs and provenance are collected at
`local_results/SynTangle_rhynchospora_comparisons.zip`; the latest run pointer is
`local_results/latest_rhynchospora_visual_run.txt`. A comparison failure can leave
an input-only PDF, which the collector explicitly labels incomplete. Input import
has lightweight tests for ambiguous IDs, low/missing quality, nonfinite values,
and out-of-range coordinates; real transfers and plots run only on Pegasus.
