# Conditional layout companion

Python API: `from syntangle.conditional_layout import optimize_conditional`.
Pass a validated fixture, complete starting LayoutState, and contiguous blocks of species IDs partitioning the existing displayed row order. Keyword budgets: `seconds=3600`, `row_seconds=20`, `sweeps=3`; callbacks: `progress(state)` and `event(record)`.

This additive heuristic addresses cases where constructing the full pair-factor MILP is expensive. It never constructs JointModel. It uses the canonical exact inversion scorer on original adjacent species pairs. A row move scores both original neighbors, including fixed outside rows; unchanged pairs are additive constants. No sampling, thinning, coordinate alteration or reassignment is performed. All original chromosomes, including unlinked ones, are retained.

The first stage uses a two-state chain dynamic program to choose whether each contiguous species block should keep its current chromosome layout or mirror it. Mirroring reverses chromosome orders within each disconnected component and flips entire chromosomes, preserving the common component order. Unary costs include every internal pair and pairwise costs include every boundary. Its optimality statement applies only to this binary block-mirror domain. This stage is skipped when explicit hard orientation constraints are present; no unsupported independence assumption is made.

The conditional search then visits genus boundary rows, followed by alternating full row sweeps. It tests whole-chromosome relocations within components, legal orientation-group flips and joint row mirrors. Later sweeps also test joint mirrors of overlapping three/five-row windows. Only strictly improving exact-scored moves are accepted, and fixed outside rows remain in boundary scoring. This does not exhaust every multi-row permutation. Neighborhoods and time limits are heuristic choices, never permanent forced-state or optimality claims.

The result includes complete reconstructed state, canonical input/output pair counts, move-scope history, evaluation count, elapsed time and explicit `global_lower_bound=None`. The final state is validated against the hard orientation basis and checked for unchanged link count and a non-increasing full objective. It is labelled a conditional heuristic, not a global certificate. It does not resume or reinterpret any previous exact-search exclusion ledger. Existing solvers remain unchanged.

Use as a separately scored companion, with frozen input hashes and provenance, while independent global searches continue. Tests compare affected-pair score deltas to the full canonical objective and block-DP costs to exhaustive block-mirror enumeration; verify evidence preservation, progress events and absence of global model construction.
