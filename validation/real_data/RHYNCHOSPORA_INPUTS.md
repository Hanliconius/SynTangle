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
