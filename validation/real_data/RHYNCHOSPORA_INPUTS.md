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
