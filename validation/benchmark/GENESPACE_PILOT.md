# GENESPACE layout pilot

Submit from Pegasus with:

```bash
bash validation/benchmark/submit_genespace_pilot.sh
```

All generation, native GENESPACE calls, optimization and collection execute
under Slurm. The preparation job generates ten deterministic cases. The array
runs one case per task, with at most five concurrent tasks on `nano`. Collection
uses `afterany`, so missing/failed/timed-out cases remain explicit.

The pilot contains six small independent-descendant cases (three species, six
ancestral chromosomes, 0/1/2 events per extant lineage, two biological seeds),
three typical-scale cases (four species, 16 ancestral chromosomes, 0/1/2
events), and one conserved Lepidoptera-like case (eight species, 31 ancestral
chromosomes, zero events). Every case has randomized chromosome order and
whole-chromosome presentation reversals. Structural events use the existing
simulator's explicit nested event plans. This small pilot is not a factorial
study or evidence of general superiority; replicated biological tests follow
only if it finds a practically useful difference.

## What is being compared

- Public input presentation.
- The installed GENESPACE 1.3.1 `pull_synChrOrd()` routine nested inside
  `riparian_engine()`. The runner evaluates the exact function definition
  extracted from the installed package; it does not approximate the method.
- The same native chromosome orders with three-start greedy whole-chromosome
  flip assistance, holding those orders fixed. This assistance is our benchmark
  code, not a claim about GENESPACE's automatic behavior or an exact flip optimum.
- Current Syntangle with its ordinary legal whole-chromosome search.

GENESPACE is evaluated using **every extant species as reference**, with
`syntenyWeight` 1 and 0.5. Every reference variant and its scores are saved.
Headline results use the best score across these variants and charge the total
ordering/flip search time across all variants, not only the winning trial.

All methods receive the same public homology IDs, coordinates, lengths and
species row order. The reference uses its supplied input chromosome order and
orientation. Query chromosome-name ordering is represented by input display
rank, avoiding informational chromosome names as a hidden ancestral shortcut.
The public bed reflects initial whole-chromosome presentation reversals before
native ordering. Physical coordinates are unchanged in the canonical evidence.
Native ordering uses coordinate rank, as its own routine does; all final states
are scored from the original physical coordinates using the common scorer.

Syntangle is free to reorder every chromosome row; GENESPACE fixes the selected
reference row's order. This is an explicit algorithmic difference, not an
equal-fixed-reference search-space comparison. Testing every reference reduces
arbitrary reference bias but does not remove that difference. We compare the
actual layout strategies, not GENESPACE orthology/synteny discovery, block
filtering or its full native renderer.

The runner never reads hidden ancestry/native-display/tangle logs. Chromosome
omissions or duplicate output chromosomes are errors. These pilot bundles have
unique homology per species; results do not establish behavior on duplicated,
ambiguous or incomplete evidence. Such cases need a separate benchmark.

## Scoring and artifacts

The common objective is the existing unweighted anchor crossing count between
adjacent species rows. Figures use the same SVG renderer with physical
within-chromosome positions for all methods. JSON files retain every final
chromosome permutation/sign assignment. GENESPACE's installed function text and
package versions are saved alongside the results.

Method timing excludes plotting. Native ordering elapsed time is recorded
inside R; R launch/import/adapter time is recorded separately. Flip assistance
and Syntangle are timed inside Python. This makes kernel timing explicit; it is
not a full end-to-end software benchmark, and tiny R timings may have limited
resolution. Each array task has a 25-minute external limit and a 30-minute Slurm
allocation. Search/node caps are not wall-clock limits. Baselines and partial
figures are checkpointed before the solver starts. Incomplete tasks remain
incomplete, never inferred as wins or proven solutions. A native baseline
beating a claimed proven Syntangle optimum fails the task as a regression.

After collection, inspect:

- `local_results/genespace_pilot/comparison_summary.md`
- `local_results/genespace_pilot/comparison_results.tsv` (scores, timings,
  status, selected reference/weight, public evidence fingerprint)
- `local_results/genespace_pilot/comparison_index.html` (linked vector gallery)
- `local_results/genespace_pilot/comparison/CASE/reference_variants.tsv`
- `logs/gs_compare.ARRAYID_TASKID.{out,err}`

Retain the entire output directory when copying the gallery: its HTML links
refer to per-case reports. Rerunning submission repeats the deterministic pilot
and replaces that pilot's result files; save prior runs before rerunning if
needed. Bulk generated results should not be committed.
