# Independent exact optimality certification experiment

This checkout is an experimental fork of SynTangle at
`2683d5cf7b802585efc2c3c1b29c6d3357a6380d`. Nothing here changes the
optimizer, plotting code, rules, or `pegasus-genespace-pilot` branch.
Python 3.10+; no Python dependencies. No images or large benchmarks are generated.

## What is implemented

* Independent validation and strict pairwise geometric scoring, using rational
  coordinates parsed directly from original JSON. No optimizer imports.
* Exact zero proof: a legal candidate with zero crossings attains the nonnegative
  lower bound.
* Exact exhaustive proof for small problems, including positive optima. Every
  permutation and every orientation satisfying hard equations is scored.
* Replay verification of exhaustive/zero certificates and SHA256 binding to input,
  candidate and checker. A saved `status` is not accepted as proof.
* Full-space integer LP export for the proposition `C <= U-1`. The exporter
  independently builds ordering, orientation and crossing constraints.
* A CPU Slurm submission helper. Tests are lightweight checks of the checker.
* [Structural lower-bound theorems](STRUCTURAL_THEOREM.md), exact cycle witnesses,
  and a constructive graph class. `structural.py` can certify positive optima
  without enumeration when a candidate attains its independently computed bound.

## What is not yet implemented

Large-case **end-to-end MILP certificate verification** is not complete. The SCIP
command file enables exact mode and proof logging, but this experiment does not
yet parse a VIPR certificate, compare its original problem to our exported model,
and verify that its conclusion proves infeasibility of this exact cutoff.
Therefore an LP export, SCIP status or successful external checker invocation
alone is never promoted to `EXACT_EXHAUSTIVE` or `EXACT_ZERO`.

SCIP 10 exact mode requires an appropriately built solver (GMP, MPFR, Boost
multiprecision and an exact LP solver). Enable exact mode **before reading** the
model. Default cutting-plane proof logs may need completion with `viprcomp`.
See [SCIP exact mode](https://www.scipopt.org/doc-10.0.0/html/EXACT.php)
and [VIPR](https://github.com/scipopt/vipr). We do not assume Pegasus has these.

This is computational proof with a reviewable checker, not a theorem-prover
formalization of the checker or its compiler/runtime. Independent replay protects
against trusting a solver's floating-point claim; software correctness remains
part of the trust base. Tests support correctness but are not its mathematical proof.

## Proof statement and scope

Let L be all legal layouts of the supplied evidence: fixed species row order,
whole-chromosome permutations, whole-chromosome flips, and supplied XOR orientation
equations. For a layout l, C(l) counts pairs of homology anchors whose **strict**
left-to-right order differs between adjacent species rows. Tied anchors do not
cross. Each homology must have at most one occurrence per species. Every supplied
chromosome and scored link is retained; no filtering occurs.

The claim is `C(candidate) = min_{l in L} C(l)`. It does not establish minimum
ribbon area, a minimum weighted objective, or minimum evolutionary rearrangements.
Unsupported extra constraint or block fields fail closed. Original decimal
coordinates are rational here: near-ties can differ from the current optimizer's
floating-point score. Such a discrepancy must be investigated before applying a
certificate to an optimizer report.

For zero: all C(l) are nonnegative integers, and a legal candidate attains zero.
For exhaustive positive proof: the iterator covers the Cartesian product of each
species' full permutations and each chromosome's two orientations, rejecting only
assignments that violate explicit hard orientation equations. If none has C<U,
U is optimal. A space or time limit always yields UNRESOLVED.

## Integer model equivalence

For every indexed chromosome pair a<b in each row, x_ab=1 means a is after b.
For each triple a<b<c, impose `0 <= x_ab+x_bc-x_ac <= 1`. These inequalities
exclude exactly the two cyclic tournaments on a triple. A tournament without a
directed triangle is transitive and hence defines a unique chromosome permutation.
All pairs are represented, including those absent from the objective.

For every chromosome use f_a=1 for reversed orientation. A hard XOR equation is
`f_a-f_b=0` for parity zero, or `f_a+f_b=1` for parity one.

For a link pair, each endpoint's strict precedence predicate is either x_ab,
1-x_ab, f_a, or 1-f_a. If endpoints tie on a chromosome, that pair contributes
zero for all layouts. Otherwise the crossing indicator is XOR of the two
predicates. Its binary truth table gives an integer polynomial
`t00 + (t10-t00)u + (t01-t00)v + (t11-t10-t01+t00)uv`.
Repeated variables are collapsed by evaluating their shared truth table.

Replace each product uv by z with `0<=z<=1`, `z<=u`, `z<=v`, and
`z>=u+v-1`. For binary u,v these force z=uv exactly. Summing crossing indicators
gives exactly C for every legal layout. Conversely, every feasible binary
assignment defines a legal permutation/orientation and these same crossings.
The exported cutoff subtracts the constant term exactly and uses integer U-1.
Thus infeasibility of that model plus a validated candidate is sufficient for
optimality. Exporting it does not establish infeasibility.

## Hypergraphs

The objective/constraint factor hypergraph connects decision variables sharing
factors or hard constraints. It is useful for exact decomposition, width measures,
elimination and reductions in the optimizer. It is **not** by itself the final
optimality certificate. This checker intentionally covers the unreduced legal
space independently; it never restarts or modifies an optimizer search. That
allows later checks of reduction correctness without assuming it at the outset.

## Run on Pegasus

Use plain saved layout JSON with `chromosome_order` and `chromosome_orientation`.
If an optimizer result nests that object, extract it explicitly first; the checker
does not guess which layout was intended.

```bash
bash certification/submit.sh check /absolute/fixture.json /absolute/layout.json
```

The helper prints the job ID and absolute result/log locations. To replay:

```bash
bash certification/submit.sh verify /absolute/fixture.json /absolute/layout.json /absolute/certificate.json
```

To prepare a larger proof obligation without making a proof claim:

```bash
bash certification/submit.sh export /absolute/fixture.json /absolute/layout.json
```

In that export directory, an exact-enabled SCIP can read `scip.commands` on stdin.
Preserve the LP, manifest, complete solver log, proof log, checker version and
checker output. Model/certificate binding and exact conclusion checks are the
next implementation milestone; no large-case certificate is certified by this release.

For a solver-independent structural check on real data, use:

```bash
bash certification/submit_structural.sh check /absolute/fixture.json /absolute/layout.json
```

This prints the result directory. It certifies optimality only when the candidate
attains the structural bound, and otherwise reports an unresolved gap. The statement
about unfinished large-case certificates above refers to the MILP/VIPR route.

### Saved planarian cases

`bash certification/submit_real_audit.sh` submits ten concurrent CPU tasks:
both planarian cases, each with five saved filtering schedules. It defaults to
the Pegasus pipeline checkout and `local_results/planarian_filter_9ol27C`.
Optional arguments specify the pipeline and filter-run directories.
The runner locates an original full-evidence fixture by the saved SHA256,
copies that original JSON without changing coordinates, validates the saved
whole-chromosome layout, and independently recomputes its strict rational score.
It does not assume the production lower bound. A scoring disagreement fails.

An attained structural bound produces `EXACT_STRUCTURAL` and is replayed in a
separate checker process. Otherwise it reports `UNRESOLVED_STRUCTURAL_GAP`.
An unresolved outcome does not invalidate a production optimum; this structural
bound can be weaker than another exact method's bound. All results, including
unresolved gaps, include hashes and are limited to the full imported fixture,
not to omitted source evidence, a weighted objective, or the original drawing.
No solver, plots, or new homology are generated by this audit.

## Isolation and publication

`exact-certification-experiment` is a separate development branch and local
checkout, not yet a separate GitHub repository. GitHub repository creation/forking
is unavailable through the connected tools. For complete repository isolation,
create an empty `Hanliconius/SynTangle-proof` repository and push this branch there.
Do not merge it into the working tool until independently reviewed and validated.
