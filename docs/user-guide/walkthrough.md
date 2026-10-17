# User guide: from synteny evidence to an improved layout

## Choose the right workflow

**Flip assistance is a helper for an existing GENESPACE plot.** GENESPACE has already discovered synteny and chosen chromosome order. `syntangle_flip_assist()` keeps that order and searches for better whole-chromosome orientations. You can use its returned inversion list in GENESPACE. It is an add-on written for SynTangle, not a native GENESPACE feature.

**Full SynTangle is a separate layout optimizer.** `syntangle_run()` searches both chromosome order and orientation using the supplied synteny evidence. GENESPACE is one possible input source; canonical tables, experimental JCVI blocks and experimental minimap2 alignments are other inputs. Full SynTangle does not run homology discovery. Its optimized order must be applied as well as its flips; an inversion list alone cannot reproduce it.

| Question | GENESPACE flip helper | Full SynTangle |
|---|---|---|
| Keep my existing chromosome order? | Yes | No: order may change within each species row |
| Change species row order? | No | No |
| Reverse whole chromosomes? | Yes | Yes |
| Move or reverse individual blocks? | Never | Never |
| Remove inconvenient links? | Never | Never |
| Search method | Strictly improving orientation changes | Adaptive neighborhoods and mixed-integer optimization |
| Optimality information | Heuristic result, no certificate | Feasible score, lower bound and remaining gap |
| Draw with GENESPACE afterwards? | Reuse original arguments with the replacement inversion list | Use SynTangle's PDF or apply the complete layout in a custom renderer |

Both use the same Python engine from R or the command line. Start with the [installation instructions](../../README.md#install), then use the [GENESPACE walkthrough](genespace.md) for your existing R objects or the [table format](inputs.md) for other evidence.

## What happens to a chromosome?

Think of each chromosome as a rigid strip with syntenic intervals marked on it. You may move the entire strip or turn it around. You may not cut it, exchange two intervals, or reverse one interval independently.

For length L, a whole-chromosome flip maps every position x to L-x. For example, a directed interval from 20 to 40 on a chromosome of length 100 becomes a directed interval from 80 to 60. The descending endpoints matter: they retain the interval's direction. Sorting those endpoints independently would conceal orientation information and is not an allowed operation.

The input tables remain unchanged. SynTangle stores a separate rank and orientation for each chromosome and applies that transformation when scoring and drawing. Relative flips are listed in `flips.tsv`; the final complete arrangement is in `layout.tsv`.

## What is being minimized?

The standard objective counts **crossings between link midpoints**. For two links between the same adjacent species rows, compare their left-to-right midpoint order on the first row with their order on the second. If the orders are strictly reversed, they contribute one crossing. Equal positions are ties and contribute zero. Crossings are summed over adjacent species pairs.

A wide block and a narrow block each contribute one link to this objective. It is not a count of genes, rearrangements, overlapping ribbon pixels, or evolutionary events. Changing the definition of the input blocks changes the problem being solved.

Moving a chromosome changes the positions of every link attached to it. Flipping it reverses the order of those positions. These operations can resolve many crossings at once while retaining the supplied evidence. They can also introduce crossings elsewhere: an interior species row connects to both neighbors. SynTangle scores the total, so one species pair can become worse while the overall objective improves. Inspect `pair_stats.tsv` to see where gains and losses occurred.

## How the search works

The flip helper evaluates permitted whole-chromosome orientation changes and accepts a strictly improving change, then repeats until it finds no improvement or reaches its budget. Where the input imposes hard orientation relationships, related chromosomes flip together. A stable answer is a local result: a combination of moves that requires an intermediate worsening may be missed.

Full SynTangle represents chromosome precedence and orientation as decisions in a mixed-integer model. Order-consistency constraints prevent impossible cyclic orders. Crossing terms connect those decisions to the objective. The implementation separates independent components, retains legal zero-crossing components, and allocates search time to the remaining components. For larger models, adaptive searches optimize windows of species rows with outside decisions fixed, then supply the improved layout to the global component search. The solver searches further and computes lower bounds. Conditional neighborhood results are not treated as global proofs.

The current feasible layout provides an upper bound on the minimum possible crossing count. Solver relaxations and search provide lower bounds. A time limit can leave a useful layout without resolving whether a better one exists. More time is not a promise of a better result in independently restarted runs; compare scores and retain the best validated layout.

## Why a low crossing score can still look tangled

A ribbon can twist between its own directed endpoints even when its midpoint crosses no other midpoint. This is why standard flip assistance and full midpoint optimization may leave visibly twisted ribbons.

Full SynTangle optionally adds ribbon refinement:

- **`crossing-cap`** searches for fewer twists while keeping crossings no higher than the layout reached by the initial crossing search. This is a cap at the current score, not necessarily at a global optimum.
- **`combined`** minimizes crossings plus `weight` times a twist metric. Some extra crossings may be accepted in exchange for fewer or narrower twists.

The `count` twist metric gives each twisted ribbon weight one. The `width` metric assigns each ribbon the smaller of its two interval widths, divided by the mean of these smaller widths over all supplied links. It emphasizes broad ribbons; its units follow the chosen input coordinates (base pairs or gene ranks). It is not a measurement of rendered pixel area.

```r
fit <- syntangle_run(
    rp, useOrder = FALSE, output = "results/full_refined",
    seconds = 600, ribbonObjective = "crossing-cap",
    twistMetric = "width", ribbonSeconds = 300
)
fit$stats
fit$pairStats
fit$pdf
```

`useOrder` must match the original GENESPACE display. See [existing inversions and replotting](genespace.md) before applying returned flips. Ribbon refinement is optional and is not part of the standard flip helper.

## Read the statistics

These are optimization diagnostics, not hypothesis tests: there are no p-values or confidence intervals here.

| Quantity | Interpretation |
|---|---|
| Input crossings | Score of the supplied starting layout |
| Output crossings U | Score of the final feasible layout; smaller is better for the crossing objective |
| Lower bound L | Solver-reported bound on the minimum for the same evidence and allowed moves |
| Gap U-L | Maximum possible further crossing reduction consistent with that bound |
| Per-pair crossings | Contributions of each adjacent pair of species to the total |
| Link/chromosome counts | Evidence and display size; these should be preserved |
| Twists | Number of ribbons whose endpoint directions oppose one another in the display |
| Width-weighted twists | Sum of normalized interval-width weights for twisted ribbons |
| Runtime and parameters | Time consumed and settings needed to interpret/reproduce the run |
| Input hashes and versions | Identify the exact input bytes and software environment |

For example, if U=300 and L=280, the encoded minimum lies between 280 and 300, subject to the validity of the solver bound. At most 20 further crossings can be removed. This does not say that 20 will be removed, or that the current layout is biologically wrong. A reduction from 500 to 300 is 40% relative to that starting layout. This percentage is not itself an optimality gap, and is undefined when the starting score is zero.

`stats.json` distinguishes three outcomes:

- `heuristic_no_certificate`: no lower bound establishing global optimality; used for flip assistance.
- `bounded_best_known`: a feasible solution with a nonzero gap.
- `solver_bound_closed`: the numerical crossing lower bound equals the feasible crossing score for the encoded model.

The public workflows set `independent_certificate=false`. Numerical bound closure is not an independently replayed exact mathematical proof. It also does not certify the biological input, a different filtering level, or the ribbon objective. After combined refinement the reported crossing gap refers to the final crossing score; the combined-objective diagnostics are in `ribbon_solver.json`.

## Save, inspect and report

Open the comparison PDF and inspect the changed chromosome signs and orders. Check the before/after totals, per-pair scores and preservation fields. A successful preservation audit checks the supplied data; it cannot establish whether upstream homology inference was correct.

Keep the entire output directory. `fixture.json`, `input_layout.tsv` and `layout.tsv` preserve the evidence and both layouts; `stats.json` records hashes, versions and settings. The [output reference](usage.md#outputs) describes every file. R also returns and saves a `syntangle_result.rds` object.

For a paper, report the source dataset and preprocessing, coordinate mode, fixed species row order, allowed moves, objective, time budget, software versions, final score, lower bound and gap. State explicitly if an example was downsampled or reconstructed rather than taken directly from published blocks. Cite the original data and any upstream synteny software separately from SynTangle.

## Real-data examples without large genome downloads

Layout optimization requires chromosome extents and directed synteny intervals, not the original FASTA assemblies. A compact example can therefore distribute the complete derived tables while linking to the original genome repository. Record the paper, exact assemblies, upstream pipeline, data license and table checksums alongside it.

If a subset is needed, document the deterministic selection rule and retain a manifest of included links. Call it a subset demonstration: its score and bound apply only to that subset. Do not remove links from a run during optimization or represent a subset as a complete reproduction of the paper figure. The current bundled [quickstart](../../examples/quickstart/README.md) is a synthetic installation example, not a published dataset.
