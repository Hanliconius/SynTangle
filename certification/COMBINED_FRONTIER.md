# Selective root cuts followed by frontier proof search

The experiment stays on the certification branch. Floating LP solutions select
violated valid four-cycle inequalities; rational multipliers provide the proof.
A root witness is retained as additional cuts are added (zero multipliers for
new rows). Frontier search starts with that exact root witness, rather than a
zero witness. Each child inherits the parent's witness. A weaker LP proposal
cannot replace a stronger exact bound. Replay verifies the branching cover,
selected valid rows, exact root witness and all node bounds. The final lower
bound must be at least the retained root bound.

The single-core Pegasus benchmark runs five variants sequentially on each of
two existing fixtures: all-cuts root, all-cuts frontier, selective root, selective
frontier with caps 101 and 501. Each gets up to 120 seconds of search; reported wall time includes preparation,
and replay time is recorded separately. No LP status alone
certifies infeasibility. The strongest independently replayed result for each
fixture is retained in strongest_verified_bounds.json, so an inferior variant
cannot replace a better certificate in the report. Graph-derived cuts do not
exclude any legal chromosome arrangement.

Submit certification/submit_combined_benchmark.sh from an isolated checkout of
exact-certification-experiment. One node and one CPU core, 8 GB, nano, 25 minutes.
The production solver is not modified. Replays hash all relevant checker files;
regenerate certificates when checker code changes.
