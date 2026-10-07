# Final validation and transition to real data

The current solver remains experimental. These jobs are final checks for the
current development round, not a declaration that the paper is fully validated.
No solver benchmark or rendering is run on the login node or by the assistant.

## Submit

From the Pegasus repository:

```bash
bash validation/benchmark/submit_final_validation.sh
```

This snapshots code with SHA256 hashes, sets up a pinned Python environment,
submits three correctness groups and six independent data-preparation tasks.
Neither array has a concurrency throttle. Existing unittests include solver,
monotone handoff, bound, mirror, deadline, saved-layout and refinement checks.
Group 1 independently enumerates legal whole-chromosome permutations/flips for
20 seeded graphs (two graph shapes), including hard orientation equations, and
compares direct, mirror and selective hybrid solves with the brute-force optimum.
It also checks checkpoint monotonicity, evidence fingerprints and zero-deadline
bounds. Group 0 includes the new tests too; group 1 isolates their runtime/failure.

```bash
cat logs/st_final_checks.*.out
tail -n 8 logs/st_final_checks.*.err
cat logs/st_real_prep.*.{out,err}
```

Only after these outputs pass, launch real-data comparisons:

```bash
bash validation/real_data/submit_block_comparisons.sh "$(cat local_results/latest_final_run.txt)"
```

All twelve tasks run without a throttle: six fixed datasets x 150/600 seconds.
Each task saves native orders, all ordering/flip variants, baselines, solver
checkpoints, graph metrics, bounds, raw candidate and guarded returned result.
The GENESPACE environment must contain the already-used 1.3.1 version and
`data.table`. Override `SCT_GENESPACE_ENV` if needed. The helper intentionally
fails if the installed version/function differs.

```bash
cat logs/st_real_compare.*.{out,err}
```

## Acceptance and limitations

- All exhaustive outcomes must equal independent brute force; every lower bound
  must be admissible, including interrupted solves.
- Imported source bytes must match pinned Git blob SHAs; both directions must
  describe the same multiset of geometric links. Otherwise import fails.
  Conflicting directional orientation annotations are preserved and flagged;
  neither is used as a hard whole-chromosome parity constraint.
- Evidence fingerprints must remain unchanged and returned states must be legal.
- The never-worse guard is a product guarantee, NOT a raw solver win. Report raw
  outcomes, ties and regressions separately.
- A completed process with unequal bounds is bounded, not proven optimal.
- A scheduler timeout or missing result is incomplete, not a loss or win.
- Repeatable numerical solver certificates are not independent rational proofs.
- Existing handoff tests cover the residual/branch path and hybrid matrix reuse;
  they do NOT prove every conceivable reduction transfers between all backends.
  No previous branch frontier is claimed to be imported into a fresh global solve.
- The current block comparison moves all chromosome rows. A matched fixed-
  reference experiment still needs an explicitly constrained solver mode before
  making claims about that setting. Do not infer it from baseline row behavior.

Before a paper: use multiple independent studies, benchmark gene-anchor and
block-level objectives separately, repeat timing on comparable nodes, and audit
exact published layouts. No solver default changes are authorized by a passing
suite alone.

## First final-check run and correction

Run final_NFyQeq (setup 138406721, checks 138406722, prep 138406723) failed the
new exhaustive test because a legal hard orientation equation crossed homology
components. Solver decision partitions now merge homology components linked by
hard equations. Layout, local/auto, branch, layer DP and residual paths share this
partition helper. Biological incidence graphs and their metrics are unchanged.
New regression tests compare a hard-coupled, homology-disconnected fixture against
independent brute force. These are pending Pegasus verification, not claimed passes.

The third check group also lacked an already-committed test module in the restored
Pegasus worktree. Rerun from a complete isolated checkout, not a partial set of
restored files. The launcher optionally copies a previous real_data preparation
folder via SCT_REUSE_PREP; all pinned hashes/import checks are repeated, but
existing source files are not downloaded again. Original failed-run files remain
unchanged.
