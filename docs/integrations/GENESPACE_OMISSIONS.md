# Why native ordering omitted chromosomes

The Rhodotorula comparison calls the nested `pull_synChrOrd` ordering function,
not the complete GENESPACE discovery/plotting workflow. The upstream implementation
in `jtlovell/GENESPACE/R/plot_riparian.R` joins reference and nonreference anchors
by orthogroup ID, computes chromosome order only for those matched rows, and uses
an inner merge of chromosome lengths with that order table. Nonreference
chromosomes with no qualifying shared reference anchor do not survive that merge.

This is not evidence of GENESPACE corrupting or rearranging a chromosome. Our
adapter expected a total order from a function that can return only a partial one.
Our wrapper's former `stop()` converted that partial result into a misleading fatal
message. Valid partial results now exit normally with `NATIVE_STATUS=PARTIAL`.

An independent check of the archived input bed/clens files against the recorded
missing chromosome IDs found exact set equality for every reference below:

| Reference | Predicted no-reference-anchor chromosomes | Reported omissions |
|---|---:|---:|
| JCM10020 | 4 | 4 |
| NBRC10512 | 4 | 4 |
| JCM10049 | 4 | 4 |
| JCM10021 | 5 | 5 |
| NBRC10513 | 5 | 5 |
| JCM10296 | 3 | 3 |
| JCM10295 | 4 | 4 |
| JCM8202 | 17 | 17 |
| JCM10449 | 1 | 1 |
| JCM10450 | 1 | 1 |
| JCM10213 | 23 | 23 |

Evidence: `validation/real_data/saved_results/SynTangle_simple_reconstructed.zip`
at commit `6f80eb5`, Rhodotorula 600-second native_input and native_variant_audit.
Both tested synteny weights have the same missing sets. No transitive orthology
or fabricated anchors are introduced to suppress these omissions.

Completion retains native relative order and appends missing chromosomes in
original input order. It is explicitly a SynTangle adapter policy, not an output
of the untouched native function. The two supplemented baseline panels and all
methods are rescored on the unchanged adjacent-link fixture. Future end-to-end
GENESPACE comparisons must be reported separately.

`submit_rhodotorula_adapter_check.sh` verifies the installed wrapper without
rerunning alignments, optimization or PDF generation. It archives the installed
ordering function and audits the resulting complete layouts.
