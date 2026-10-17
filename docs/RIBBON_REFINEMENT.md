# Ribbon-aware refinement

`refine-ribbons` refines a complete supplied chromosome layout while retaining all
fixed block coordinates, strands and links. Only whole-chromosome permutations
and flips are allowed. Install with `pip install '.[optimization]'`.

```bash
syntangle refine-ribbons fixture.json --layout layout.tsv --output refined --seconds 600
```

The default minimizes width-weighted self-twists subject to no increase in the
starting midpoint crossing count. This is a crossing cap, not a claim of
lexicographic global optimality. Run the standard crossing optimizer first when
crossing minimization is the primary goal.

To allow a crossing/twist trade-off:

```bash
syntangle refine-ribbons fixture.json --layout layout.tsv --output balanced --objective combined --weight 1 --seconds 600
```

Combined minimizes `C + weight*T`. Select `--twist-metric count` for unweighted
self-twists; the default `width` uses the smaller endpoint block span divided by
the mean such span over all links. This is a visual proxy, not a biological cost
or an exact overlapping-area metric. A self-twist means a ribbon's two endpoint
orders disagree after whole-chromosome orientations are applied.

Input TSV columns are `species`, `chromosome`, `rank`, `orientation` (1 or -1).
Every fixture chromosome must appear exactly once; orientation constraints apply.
The same format is accepted by the GENESPACE refinement interface.

Outputs are `layout.tsv` and `audit.json`: both crossing and twist metrics, input
hashes, allowed moves, solver status, modified-objective numerical dual bound and
full final state. If no better feasible solution is found, the starting layout
is retained. Numerical solver bounds are not independent exact certificates.
Legacy optimization commands and their midpoint objective remain unchanged.

The Rhodotorula evaluation retained 288 midpoint crossings while reducing twists
from 211 to 177 with the width cap; the combined width trial reached 308 crossings
and 154 twists. These are case-specific feasible results, not universal settings
or new optimality claims.
