# Paired local-search scoring benchmark

Run `submit_scoring_benchmark.sh` from a separate checkout while the original
stress array continues. Input bundles are read from the original stress root;
outputs and code are isolated. Each of nine jobs runs the old and cached local
search sequentially, one restart and five improving steps with the same seed,
120 seconds per method, and a six-minute Slurm limit. All nine may run at once.

This measures the diagnosed local-search bottleneck, not total exact-solver
runtime. Cached timings include term construction and validated progress writes.
Both completed methods must produce the same complete layout and diagnostics
(including candidate and subset-DP evaluation counts) before a speedup is claimed.
A timeout is reported explicitly; no incomplete ratio or optimality claim is made.
Original GENESPACE comparison results remain in their original directories.

The scorer groups pairs of homology links by their chromosome-order/orientation
decisions and tabulates their exact four-entry costs. Candidate scores update only
terms whose decisions change. Same-position ties and reflected floating-point
coordinates use the original anchor-key operations. Candidate moves, biological
constraints, randomization, and tie-breaking remain the same. The final selected
score is checked with the canonical scorer.

The automatic branch-and-bound fallback now forwards `local_max_improving_steps`
to its preliminary local search (previously ignored). The default remains 10000.
An optional callback exports strictly improving local-search incumbents; comparison
runs save these atomically to `incumbent.json`. These are heuristic layouts, not
completed exact results, and they are not classified as completed comparisons.
Branch-and-bound improvements after the preliminary search are not yet streamed.
