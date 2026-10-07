# Mirror equivalence and selective neighborhood experiment

The efficiency ablation completed all 252 tasks. All 36 largest-case runs at
600 seconds reported C=lower=8377, including public starts. Adaptive neighborhoods
improved public-start layouts at 150 seconds by about 69–82% relative to plain
weighted no-hint MILP. Timing rankings varied; six strict-improvement audits
verified smaller cases, while three largest audits were unresolved at 120 seconds.
See the [archived evidence](results/2026-10-06-efficiency/summary.md).

This next experiment keeps the default production solver unchanged. It adds one
proved equivalence reduction and one explicitly heuristic scheduling policy.

## Exact mirrored-state reduction

For each component the joint objective is

    E(x) = c + sum_i l_i*x_i + sum_(i,j) q_ij*x_i*x_j.

The implementation verifies, with integer arithmetic and no tolerance,

    2*l_i + sum_j q_ij = 0 for every primary decision i.

This identity establishes E(1-x)=E(x). The feasible primary domain is also closed
under full complement: row precedence complement reverses a complete total
order, preserving all triangle transitivity constraints; free GF(2) group
complement preserves the hard orientation constraints. Product variables are
reconstructed from complemented primary bits. This proves a two-state symmetry
orbit and permits retaining only one representative.

Only an unrestricted direct global component solve applies this reduction.
It anchors the primary bit with the largest sum of incident absolute quadratic
coefficients (deterministic index tie-break) to the current incumbent's value so
that incumbent and
its MIP start remain feasible. The column fixing persists throughout the single
global solver run; there is no subsequent backend transition that restores it.
Every result records the exact check, fixing and reconstruction rule. A failed
identity or noninteger coefficient disables pruning. Existing conditional
neighborhood calls are unchanged; their fixings never become global exclusions.

This removes a factor of two from the represented primary search space. It
**does not promise a factor-of-two speedup**: presolve may already recognize
symmetry, or LP/cut work may dominate. Both possibilities require measurement.
A future objective or hard-domain change must revisit the certificate before
using this reduction. Full-layout reconstruction remains legal and reversible.

## Selective neighborhood policy

`mirror_selective` tests a compound policy:

- Skip neighborhoods for components with at most 512 primary decisions.
- Spend at most 15% of the component allocation and at most 30 seconds on them.
- Use adaptive three-to-five-species neighborhoods for allowances below 300 s;
  use three-species neighborhoods for longer allowances.
- Visit each distinct adjacent-species window once per sweep, avoiding repeated
  clamped windows at row ends. A new sweep can revisit it after improvements.

These are exploratory thresholds, not established biological scaling laws or
certified pruning. Skipping a neighborhood leaves every unresolved choice
available to global search. An unchanged sweep is still heuristic stagnation,
not proof that its conditional region cannot improve. Stopping this work retains
more time for the one global solve.

## Controls and design

Seven variants: direct no-hint, three-species and adaptive controls; each with
mirror equivalence; and mirror_selective. The six paired variants isolate the
mirror fixing. The selective variant bundles four policy choices and cannot
attribute a gain to any one of those choices without a later ablation.

All variants use the same direct HiGHS version, one CPU, 16G, input fingerprints
and frozen public/pre-hybrid starts from the earlier efficiency run. Reference
starts are copied and re-scored, not updated to the now-known optimum. The
launcher finds the most recent 252-task efficiency directory; override with
`SCT_REFERENCE_RUN=/absolute/path/to/efficiency_RUN` when necessary.

All nine presentations x two starts x seven variants at 150 s once = 126 tasks.
The three largest x two starts x seven variants at 600 s twice = 84 tasks.
Total 210, submitted without a throttle. Longer repeats address timing variability
seen in the previous ablation. Incumbents, memory, stage costs, nodes and mirror
fixings are saved. The workers verify input fingerprints before optimization.

A separate nine-task strict-improvement audit keeps the unanchored formulation,
now allowing 600 seconds per positive-score component. It re-scores historical
reported optima and remains explicitly unresolved if that allowance is exhausted.
Audit uses the same numerical backend, not an independent rational proof checker.
No outcome or new default is inferred before these runs complete.

## Run on Pegasus

    bash validation/benchmark/submit_policy_benchmark.sh

Read the final collector and audit logs:

    cat logs/st_policy_collect.*.{out,err}
    cat logs/st_policy_audit.*.{out,err}

Source is snapshotted at submission. Job IDs, tested commit, frozen starts and
input hashes are retained in the unique `local_results/policy_*` directory.
