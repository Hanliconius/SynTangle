# Selective cycle proof experiment

The previous single-core Pegasus benchmark gave planarians U=L=7901 at one root with four-cycle inequalities. For flatworms Fig. 4f, baseline L=43879 increased to 53345 using the root alone (5.043s), and to 55339 with best-bound search (101 nodes, 50.914s); U remained 60294. Every certificate replayed. The hard case is unresolved.

This next experiment compares all 4000 root cuts against iterative batches of 16 or 64 violated cuts from the same graph-selected pool. It uses one core and one search node per variant, with a 60s separation budget. LP calls and wall time, including replay time, are reported separately. Fewer rows do not guarantee faster solves: repeated LP construction may offset savings. Measurements are pending.

Cut selection uses floating LP solutions only to rank useful valid inequalities. The verifier regenerates each chosen inequality from the original model and independently replays rational multipliers. The feasible integer space is unchanged. Every previously certified root witness is retained, extended with zeros when new rows are added, and replaced only by an equal or stronger exact bound. These are fresh comparisons, not continuation of previously saved trees.

`bash certification/submit_selective_benchmark.sh` reads the existing two saved real fixtures, runs three variants sequentially per fixture, and replays all six certificates. No production optimization code changes.
