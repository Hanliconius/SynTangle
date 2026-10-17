# Methods and guarantees

A layout assigns a permutation and an orientation sign to each species' chromosomes. Species row order and all within-chromosome evidence are fixed. A negative sign maps every coordinate on a chromosome of length L to L-x. Internal rearrangements and independent ribbon endpoint sorting are prohibited.

For each pair of adjacent species, a link is scored at the midpoint of its fixed interval on each chromosome. Two links cross when their midpoint orders are strictly reversed across those rows. Ties do not count as crossings. The total objective is the sum of these strict inversions across adjacent species pairs. Rendering a thick ribbon introduces visual features that this midpoint objective does not count.

Standard flip assistance repeatedly evaluates legal orientation changes and accepts a strictly improving change, keeping the chromosome permutation fixed. This is a heuristic: a locally stable result need not be globally optimal.

Full SynTangle combines heuristic improvement and mixed-integer optimization over chromosome ordering and orientation. It retains a feasible incumbent and a crossing lower bound. A feasible score U and valid bound L give L <= optimum <= U. Numerical equality closes the solver-reported model gap; an independently replayable exact certificate requires additional verification and is not emitted by these public workflows.

## Ribbon objective

A ribbon is twisted when its directed endpoints have opposite orientations on the displayed chromosomes. The count metric weights every ribbon equally. The width metric uses the smaller of its two interval widths, divided by the mean such width across the supplied links.

Crossing-cap refinement minimizes twists subject to C <= C_start, where C_start is the score after the initial crossing search. Combined refinement minimizes C + weight*T. This can exchange some midpoint crossings for fewer or narrower twisted ribbons. All moves remain whole-chromosome moves; all evidence is retained.

## Scope and reproducibility

Scores describe the exact supplied evidence, coordinate system, species order and allowed moves. They cannot be compared directly across different block definitions or filtering levels. Always report the input source, preprocessing, coordinate mode, solver/software versions, time budget, objective, final score and bound. Preserve the emitted hashes, fixture, layout and diagnostics.

These are layout methods, not estimates of evolutionary event counts. Synteny evidence reconstructed from genomes is a reanalysis and should not be described as an exact reproduction of an author's discovery pipeline without matching the source data and parameters.
