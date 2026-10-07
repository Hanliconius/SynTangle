# Intrinsic and visual tangledness: testable contracts

The scientific target is to demonstrate that a chromosome-synteny presentation
can contain many avoidable crossings while its underlying legal layout problem
has few or none. These metrics count unweighted homology-link crossings between
adjacent species layers. They are not evolutionary rearrangement counts or a
measurement of how readers perceive a figure. The current solver concerns
chromosome-synteny layouts, not phylogenetic tanglegrams.

For the same biological evidence, fixed species layer order and hard constraints,
let D be the displayed crossing count, C* the minimum legal count, U the best
known feasible count, and L a valid global lower bound. A legal display is also
an upper bound: replace U by min(D, U).

- Visual tangledness is D.
- Intrinsic tangledness lies in [L, U]; report its exact value only if L = U.
- Excess presentation tangledness lies in [D - U, D - L].
- At least D - U displayed crossings are demonstrably avoidable, even without a proof of optimality.
- The best retained layout's excess lies in [0, U - L].

A feasible zero-crossing layout proves C* = 0 by nonnegativity. A timeout with
U > L establishes neither optimality nor impossibility of solving the graph.
The helper in `src/syntangle/tangledness.py` performs interval arithmetic; it
does not manufacture a certificate or establish that a supplied lower bound
is globally valid. Experimental hybrid results now include these metrics.
Legacy `excess_layout_crossings_initial` fields still mean crossings removed;
when a solve is unresolved, interpret them as demonstrated avoidability, not
exact excess.

## Reference-layout refinement

`src/syntangle/refinement.py` validates complete reference/candidate layouts,
including hard orientation equations, then scores both against the same fixture.
A strict improvement replaces the reference; a tie, worse candidate or absent
candidate retains it. The raw candidate score and selection outcome remain
visible. Retaining the baseline is not counted as a win.

`refine_reference_layout` starts one fresh experimental hybrid run from a copy
of the reference. It checks unchanged input evidence and canonically rescores
the returned layout. This supplies an incumbent, not an imported search frontier
or a certificate of carrying reductions across independent runs. It is an
opt-in adapter; production defaults and existing independent GENESPACE pilot
results are unchanged. It does not run native GENESPACE. All chromosome rows
remain movable, including the row used as a reference by GENESPACE; a comparison
with a fixed reference row requires a correspondingly constrained solver model.

## Tests and remaining evidence

`tests/test_tangledness_refinement.py` checks a synthetic 56-crossing display
with a proven zero optimum, intervals against an exhaustively enumerated nonzero
optimum, invalid bounds, legal-layout validation, tie/regression retention,
canonical score checks and a real hybrid refinement. Synthetic references are
explicitly identified; they are not native GENESPACE outputs.

These are software and mathematical-contract tests. The paper's claim of
consistent usefulness still needs independently generated biological instances
and real inputs, native GENESPACE and separately labelled flip-assisted baselines,
matched legal freedoms, documented total runtime budgets, and retained ties and
failures. It also needs paired presentations of identical biology to measure
avoidable drawing artifacts. Structural solvability analysis should record
reduced component size, graph width/separators, decision/factor counts and shared
ordering conflicts alongside bounds, runtime and proof status. Those experiments
are future work, not established by these tests. Pending mirror/policy experiment
results must be reviewed before adopting their defaults.

Run the contracts with:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p test_tangledness_refinement.py -q
```
