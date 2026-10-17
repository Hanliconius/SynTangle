# Opt-in ribbon-twist objectives

These experiments do not replace the midpoint-crossing objective or amend RULES.md.
All original coordinates and links remain fixed. Legal moves are whole-chromosome
order and orientation changes satisfying the same hard constraints.

For one link, let b1,b2 be its occurrence strand bits and f1,f2 the chromosome flip
bits. Its twist indicator is b1 XOR b2 XOR f1 XOR f2. This matches the relative
direction of its two directed endpoint pairs. Chromosome offsets/order do not
change this self-twist indicator. It does not count overlapping ribbon area.

Four trials start from the saved Rhodotorula 288-crossing layout:

1. Minimize twist count subject to crossings <= 288.
2. Minimize width-weighted twists subject to crossings <= 288.
3. Minimize crossings + twist count (lambda=1).
4. Minimize crossings + width-weighted twists (lambda=1).

A link's width weight is the smaller of its two genomic block spans divided by
the mean of those smaller spans over all fixture links. Thus mean link weight is
one and broad blocks matter more than narrow ones. This is a display proxy, not
an evolutionary or biological penalty. Lambda=1 is an initial trial, not calibrated.

The cap trials preserve the incumbent crossing score or improve it. They implement
a tie-break among layouts at/below the supplied cap; they are lexicographic only
when that cap is already a true minimum for the same feasible space. The archived
288 result is solver-reported closed-gap, not independently certified here.

The standalone experiment constructs one global order/orientation model, adds twist
terms and uses a 600-second, single-core SciPy/HiGHS solve. It retains the supplied
incumbent if no better feasible solution is returned. It does not inherit old
bounds/certificates as proof of the new objective. Reports include both midpoint
crossings and twist metrics, numerical solver status and the modified-objective
bound. Every experimental PDF receives the same directed endpoint checks.

`submit_ribbon_experiment.sh` submits four unthrottled tasks using cached inputs.
It writes to a new result directory and never overwrites the existing comparisons.
