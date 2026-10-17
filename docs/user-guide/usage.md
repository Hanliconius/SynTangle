# Running SynTangle and reading results

For an explanation of the workflows, allowed moves, search and statistics, begin with the [user guide](walkthrough.md). This page is the command and output reference.

`syntangle flip-assist` and `syntangle run` accept either `--chromosomes chromosomes.tsv --links links.tsv` or `--fixture fixture.json`. Both require `--output DIRECTORY`. Use `--figure` for a vector PDF. Run `syntangle --version` or each command's `--help` for the installed interface.

Flip assistance uses strictly improving whole-chromosome flips and preserves chromosome order. Full optimization uses the adaptive hybrid search. `--seconds` is the search budget (60 seconds by default for flips, 600 for full optimization); parsing, model construction and output can add overhead. `--seed` defaults to 1. Equal budgets are not a guarantee of identical results across machines or solver versions.

Optional full-search refinement:

| Option | Meaning |
|---|---|
| `--ribbon-objective none` | Standard midpoint-crossing objective (default) |
| `--ribbon-objective crossing-cap` | Minimize twists while retaining or improving the crossing count reached by the first search |
| `--ribbon-objective combined` | Minimize crossings plus `--weight` times the selected twist metric; crossings may increase |
| `--twist-metric width` | Weight twists by interval width (default) |
| `--twist-metric count` | Count twisted ribbons equally |
| `--ribbon-seconds` | Additional refinement budget, default 600 seconds |

## Outputs

| File | Contents |
|---|---|
| `layout.tsv` | Final chromosome ranks and signs |
| `input_layout.tsv` | Starting ranks and signs |
| `fixture.json` | Fixed normalized input evidence |
| `stats.json` | Before/after metrics, bounds, gap, runtime, parameters, dependency versions, hashes and preservation checks |
| `summary.tsv` | Compact statistics, also returned to R |
| `pair_stats.tsv` | Crossings for each adjacent species pair |
| `flips.tsv` | Chromosomes flipped relative to the supplied starting display |
| `solver.json` | Full crossing-search diagnostics, when applicable |
| `ribbon_solver.json` | Optional ribbon-refinement diagnostics |
| `flip_assist_comparison.pdf` or `run_comparison.pdf` | Input/output drawing with fixed link endpoints |

The R functions additionally save `syntangle_result.rds` and return an absolute inversion list suitable for GENESPACE. Existing nonempty output directories are rejected.

## Statistical interpretation

`input` and `output` in `stats.json` report midpoint `crossings`, `pair_crossings`, `link_count`, `twists`, and `width_weighted_twists`. `crossing_lower_bound` bounds the best achievable crossing count for the supplied model. `crossing_gap = output.crossings - crossing_lower_bound`. Bounds are absent for flips-only assistance.

`crossing_status` is one of:

- `heuristic_no_certificate`: an improved or unchanged layout, without a bound establishing optimality.
- `bounded_best_known`: a feasible layout and a lower bound with a remaining gap.
- `solver_bound_closed`: matching numerical solver bound and feasible crossing count for this model.

`independent_certificate` is false for these workflows. Independent proof experiments are maintained separately. Historical raw solver diagnostics may use the wording “proven optimum”; the public `stats.json` classification states the evidence level precisely. A crossing certificate does not certify a twist objective, biological correctness, or arbitrary changes to the input graph.
