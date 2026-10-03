# Experimental ancestral-state pruning

This directory is an intentionally isolated proof-of-concept. It is **not part
of the production SynTangle solver** and should not be merged into the main
pipeline unless the pruning logic is shown to preserve the intended biological
and optimization guarantees.

## Question

Can information about ancestral chromosome structure, inferred jointly from:

- the extant multispecies homology hypergraph / ordered block occurrences; and
- a rooted species tree,

provide useful lower bounds early enough to reject some residual topology
hypotheses before expensive layout search?

The experiment starts with the smallest useful version of that question.

## Representation

The biological data remain the SynTangle multispecies hypergraph represented by
the chromosome-homology incidence structure.

From the ordered block occurrences we derive a second object: binary adjacency
characters. For a pair of homology groups H1 and H2 in a species:

- 1: H1 and H2 occur exactly once and are adjacent on the same chromosome;
- 0: H1 and H2 occur exactly once but are not adjacent;
- ?: either homology group is missing or duplicated, so no binary claim is
  made.

An unordered block adjacency is used because whole-chromosome reversal does not
change whether two blocks are neighbors.

## Tree calculation

For each adjacency character, a two-state Sankoff dynamic program is run on the
rooted species tree.

For each candidate ancestral root state it reports the minimum number of
adjacency changes required.

### 1. A derived terminal rearrangement

Fixture: fixtures/terminal_fusion.json

~~~text
sp1  H1 | H2    absent
sp2  H1 | H2    absent
sp3  H1 | H2    absent
sp4  H1-H2      present

tree: ((sp1,sp2),(sp3,sp4));
~~~

Expected Sankoff costs:

~~~text
ancestral H1-H2 absent   -> 1 change
ancestral H1-H2 present  -> 2 changes
~~~

If we already possess a chromosome-valid feasible history costing one event,
any candidate root topology that requires H1-H2 to be present has an
independent-character lower bound of two events and can be rejected under that
event objective.

### 2. Genuine ancestral ambiguity

Fixture: fixtures/ambiguous_root.json

~~~text
sp1  absent
sp2  absent
sp3  present
sp4  present

tree: ((sp1,sp2),(sp3,sp4));
~~~

Expected root-conditioned costs:

~~~text
root absent   -> 1 change
root present  -> 1 change
~~~

Neither ancestral state can be removed. This is an explicit control against
turning a preferred reconstruction into a false hard constraint.

### 3. Species-tree dependence

The same extant states from the ambiguous fixture require two changes on:

~~~text
((sp1,sp3),(sp2,sp4));
~~~

rather than one on:

~~~text
((sp1,sp2),(sp3,sp4));
~~~

So the species tree is supplying real information; the procedure is not simply
taking a majority vote over extant species.

## Why this is a lower bound rather than an ancestral genome reconstruction

Adjacency characters are optimized independently in this first experiment.
The individually cheapest ancestral states need not combine into a physically
valid set of chromosomes.

That is useful rather than fatal. Independent optimization is a relaxation, so
its event count can be used as a **lower bound** on a chromosome-valid
ancestral history.

For a proposed ancestral topology T:

~~~text
LB_history(T) =
    sum of root-conditioned minimum costs
    over the adjacency states required by T
~~~

Unspecified adjacencies retain their unconditional minimum.

A topology is called prunable in this experiment only when:

~~~text
LB_history(T) > cost(best known feasible chromosome-valid history)
~~~

A singleton parsimony-preferred root state by itself is **not** treated as a
hard biological fact.

This distinction is important for eventual integration with SynTangle and its
monotone-pruning rule.

## Run the simple examples

~~~bash
python -m experiments.ancestral_pruning.ancestral_pruning_poc \
  experiments/ancestral_pruning/fixtures/terminal_fusion.json \
  --tree '((sp1,sp2),(sp3,sp4));'
~~~

and:

~~~bash
python -m experiments.ancestral_pruning.ancestral_pruning_poc \
  experiments/ancestral_pruning/fixtures/ambiguous_root.json \
  --tree '((sp1,sp2),(sp3,sp4));'
~~~

The unit tests also compare root hypotheses and the tree-dependent control:

~~~bash
python -m unittest tests.test_ancestral_pruning_experiment -v
~~~

## What this does and does not test

This first proof-of-concept tests:

- extraction of conservative adjacency characters from the existing
  hypergraph-backed fixture representation;
- rooted-tree dynamic programming;
- root-conditioned event lower bounds;
- a case where a candidate ancestral state becomes safely dominated by a known
  feasible history;
- a case where ambiguity must be retained;
- sensitivity to species-tree topology.

It does **not** yet test:

- complete chromosome-valid ancestral genome assembly;
- coupled adjacency compatibility;
- fusion versus fission direction as separate event types;
- inversions;
- branch-specific rates or asymmetric event costs;
- uncertain species trees;
- direct pruning of the production layout branch-and-bound;
- whether the historical lower bound actually improves runtime on nontrivial
  SynTangle kernels.

## Next experiment if the simple cases behave correctly

The next useful step is still small:

1. simulate 4-6 species with 6-10 chromosomes;
2. introduce one, then two, then overlapping fusion/fission events;
3. retain the known species tree and hidden ancestral genome;
4. derive adjacency lower bounds using only extant public evidence;
5. enumerate the small residual topology space both with and without the
   ancestral lower bound;
6. compare how many candidate topologies are rejected;
7. verify that the hidden true ancestral state is never incorrectly removed;
8. deliberately include non-identifiable histories and verify that alternatives
   remain.

Only after that should ancestral bounds be wired into the main residual
factor/branch-and-bound pipeline.
