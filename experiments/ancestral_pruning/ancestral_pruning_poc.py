from __future__ import annotations

"""Small proof-of-concept for tree-aware ancestral lower bounds.

This module is intentionally outside src/syntangle. It explores whether
ancestral information derived from the extant multispecies homology
representation can reduce a later layout/topology search without changing the
production solver.

The experiment treats each unordered adjacency between homologous blocks as a
binary character on a rooted species tree:

    1 = the two homology groups are adjacent on one chromosome
    0 = both occur exactly once in the species but are not adjacent
    ? = missing or duplicated evidence, so the state is unknown

A Sankoff dynamic program gives the minimum number of changes compatible with
the observations, both unconditionally and conditional on either root state.

These per-adjacency costs are an independent-character relaxation. Summing
them gives a lower bound on any chromosome-valid ancestral history. Therefore
a candidate root topology can only be pruned safely when its conditioned lower
bound is already worse than a known feasible history/incumbent under the same
event objective. A singleton best root state is not, by itself, treated as a
hard biological truth.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from syntangle import Fixture, load_fixture


Adjacency = tuple[str, str]


@dataclass(frozen=True)
class TreeNode:
    name: str | None
    children: tuple["TreeNode", ...]

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def leaf_names(self) -> tuple[str, ...]:
        if self.is_leaf:
            if self.name is None:
                raise ValueError("Leaf without a name")
            return (self.name,)
        return tuple(
            name
            for child in self.children
            for name in child.leaf_names()
        )


@dataclass(frozen=True)
class AdjacencyEvidence:
    adjacency: Adjacency
    observations: dict[str, int | None]


@dataclass(frozen=True)
class SankoffResult:
    adjacency: Adjacency
    minimum_changes: int
    root_cost_absent: int
    root_cost_present: int

    @property
    def optimal_root_states(self) -> tuple[int, ...]:
        best = min(self.root_cost_absent, self.root_cost_present)
        return tuple(
            state
            for state, cost in (
                (0, self.root_cost_absent),
                (1, self.root_cost_present),
            )
            if cost == best
        )

    def root_cost(self, state: int) -> int:
        if state == 0:
            return self.root_cost_absent
        if state == 1:
            return self.root_cost_present
        raise ValueError("Binary ancestral state must be 0 or 1")

    def excess_cost(self, state: int) -> int:
        return self.root_cost(state) - self.minimum_changes


@dataclass(frozen=True)
class TopologyBound:
    unconditional_lower_bound: int
    conditioned_lower_bound: int
    excess_lower_bound: int
    incumbent_event_cost: int | None

    @property
    def prunable(self) -> bool:
        return (
            self.incumbent_event_cost is not None
            and self.conditioned_lower_bound > self.incumbent_event_cost
        )


def canonical_adjacency(a: str, b: str) -> Adjacency:
    if a == b:
        raise ValueError("Self-adjacencies are not used in this experiment")
    return tuple(sorted((a, b)))


def parse_newick(text: str) -> TreeNode:
    """Parse the small rooted Newick subset needed by the experiment.

    Leaf/internal names are supported. Branch lengths are accepted but ignored.
    Quoted labels, comments, and Newick metadata are intentionally out of scope
    for this proof-of-concept.
    """

    source = "".join(text.split())
    if source.endswith(";"):
        source = source[:-1]
    if not source:
        raise ValueError("Empty Newick string")

    index = 0

    def parse_label() -> str | None:
        nonlocal index
        start = index
        while index < len(source) and source[index] not in ",():":
            index += 1
        label = source[start:index]
        if index < len(source) and source[index] == ":":
            index += 1
            while index < len(source) and source[index] not in ",()":
                index += 1
        return label or None

    def parse_node() -> TreeNode:
        nonlocal index
        if index >= len(source):
            raise ValueError("Unexpected end of Newick string")

        if source[index] == "(":
            index += 1
            children: list[TreeNode] = []
            while True:
                children.append(parse_node())
                if index >= len(source):
                    raise ValueError("Unclosed Newick group")
                if source[index] == ",":
                    index += 1
                    continue
                if source[index] == ")":
                    index += 1
                    break
                raise ValueError(
                    f"Unexpected Newick token {source[index]!r} at {index}"
                )
            name = parse_label()
            return TreeNode(name=name, children=tuple(children))

        name = parse_label()
        if name is None:
            raise ValueError(f"Expected leaf label at position {index}")
        return TreeNode(name=name, children=())

    root = parse_node()
    if index != len(source):
        raise ValueError(
            f"Unexpected trailing Newick text at position {index}: "
            f"{source[index:]!r}"
        )

    leaves = root.leaf_names()
    if len(leaves) != len(set(leaves)):
        raise ValueError("Species-tree leaf names must be unique")
    return root


def _species_homology_layout(
    fixture: Fixture,
) -> tuple[
    dict[str, dict[str, int]],
    dict[str, set[Adjacency]],
]:
    counts: dict[str, dict[str, int]] = {
        species: {} for species in fixture.species_ids
    }
    observed_adjacencies: dict[str, set[Adjacency]] = {
        species: set() for species in fixture.species_ids
    }

    for chromosome in fixture.chromosomes:
        species = chromosome.ref.species_id
        ordered = tuple(
            sorted(
                chromosome.blocks,
                key=lambda block: (block.start, block.end, block.occurrence_id),
            )
        )
        for block in ordered:
            counts[species][block.homology_id] = (
                counts[species].get(block.homology_id, 0) + 1
            )
        for left, right in zip(ordered, ordered[1:]):
            if left.homology_id == right.homology_id:
                continue
            observed_adjacencies[species].add(
                canonical_adjacency(
                    left.homology_id,
                    right.homology_id,
                )
            )

    return counts, observed_adjacencies


def derive_adjacency_evidence(
    fixture: Fixture,
    *,
    candidate_adjacencies: Iterable[Adjacency] | None = None,
) -> tuple[AdjacencyEvidence, ...]:
    """Derive conservative binary adjacency observations from a fixture.

    If either homology group is absent or duplicated in a species, the
    observation is marked unknown rather than coerced to absence.
    """

    counts, observed = _species_homology_layout(fixture)

    if candidate_adjacencies is None:
        candidates = {
            adjacency
            for species_adjacencies in observed.values()
            for adjacency in species_adjacencies
        }
    else:
        candidates = {
            canonical_adjacency(a, b)
            for a, b in candidate_adjacencies
        }

    evidence: list[AdjacencyEvidence] = []
    for adjacency in sorted(candidates):
        a, b = adjacency
        observations: dict[str, int | None] = {}
        for species in fixture.species_ids:
            count_a = counts[species].get(a, 0)
            count_b = counts[species].get(b, 0)
            if count_a != 1 or count_b != 1:
                observations[species] = None
            else:
                observations[species] = int(
                    adjacency in observed[species]
                )
        evidence.append(
            AdjacencyEvidence(
                adjacency=adjacency,
                observations=observations,
            )
        )
    return tuple(evidence)


def _validate_tree_species(tree: TreeNode, fixture: Fixture) -> None:
    leaves = set(tree.leaf_names())
    species = set(fixture.species_ids)
    if leaves != species:
        missing = sorted(species - leaves)
        extra = sorted(leaves - species)
        raise ValueError(
            "Species tree leaves do not match fixture species; "
            f"missing={missing}, extra={extra}"
        )


def sankoff_binary(
    tree: TreeNode,
    evidence: AdjacencyEvidence,
    *,
    change_cost: int = 1,
) -> SankoffResult:
    if change_cost < 0:
        raise ValueError("change_cost must be nonnegative")

    infinity = 10**12

    def visit(node: TreeNode) -> tuple[int, int]:
        if node.is_leaf:
            assert node.name is not None
            observation = evidence.observations.get(node.name)
            if observation is None:
                return (0, 0)
            if observation == 0:
                return (0, infinity)
            if observation == 1:
                return (infinity, 0)
            raise ValueError(
                f"Invalid binary observation for {node.name}: {observation}"
            )

        child_costs = [visit(child) for child in node.children]
        output: list[int] = []
        for parent_state in (0, 1):
            total = 0
            for absent_cost, present_cost in child_costs:
                total += min(
                    absent_cost
                    + (change_cost if parent_state != 0 else 0),
                    present_cost
                    + (change_cost if parent_state != 1 else 0),
                )
            output.append(total)
        return output[0], output[1]

    absent, present = visit(tree)
    minimum = min(absent, present)
    return SankoffResult(
        adjacency=evidence.adjacency,
        minimum_changes=minimum,
        root_cost_absent=absent,
        root_cost_present=present,
    )


def analyze_fixture(
    fixture: Fixture,
    tree: TreeNode,
) -> tuple[SankoffResult, ...]:
    _validate_tree_species(tree, fixture)
    return tuple(
        sankoff_binary(tree, evidence)
        for evidence in derive_adjacency_evidence(fixture)
    )


def topology_event_lower_bound(
    analyses: Iterable[SankoffResult],
    *,
    root_assignment: dict[Adjacency, int] | None = None,
    incumbent_event_cost: int | None = None,
) -> TopologyBound:
    """Lower-bound a proposed ancestral root topology.

    root_assignment may be partial. Unspecified adjacencies use their
    unconditional minimum. Because characters are optimized independently,
    this is a relaxation of the chromosome-valid history problem and therefore
    a lower bound rather than a complete ancestral reconstruction.
    """

    assignment = {
        canonical_adjacency(*adjacency): state
        for adjacency, state in (root_assignment or {}).items()
    }

    unconditional = 0
    conditioned = 0
    seen: set[Adjacency] = set()

    for result in analyses:
        adjacency = result.adjacency
        seen.add(adjacency)
        unconditional += result.minimum_changes
        if adjacency in assignment:
            conditioned += result.root_cost(assignment[adjacency])
        else:
            conditioned += result.minimum_changes

    unknown = sorted(set(assignment) - seen)
    if unknown:
        raise ValueError(
            "Root assignment contains adjacencies not present in analyses: "
            + ", ".join(f"{a}--{b}" for a, b in unknown)
        )

    if incumbent_event_cost is not None and incumbent_event_cost < 0:
        raise ValueError("incumbent_event_cost must be nonnegative")

    return TopologyBound(
        unconditional_lower_bound=unconditional,
        conditioned_lower_bound=conditioned,
        excess_lower_bound=conditioned - unconditional,
        incumbent_event_cost=incumbent_event_cost,
    )


def summarize(
    fixture: Fixture,
    tree: TreeNode,
) -> str:
    analyses = analyze_fixture(fixture, tree)
    lines = [
        f"fixture\t{fixture.fixture_id}",
        f"species\t{len(fixture.species_ids)}",
        f"candidate_adjacencies\t{len(analyses)}",
        "",
        "adjacency\troot0_cost\troot1_cost\tminimum\toptimal_root_states",
    ]
    for result in analyses:
        a, b = result.adjacency
        states = ",".join(map(str, result.optimal_root_states))
        lines.append(
            f"{a}--{b}\t{result.root_cost_absent}\t"
            f"{result.root_cost_present}\t{result.minimum_changes}\t{states}"
        )
    return "\n".join(lines)


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("fixture")
    parser.add_argument(
        "--tree",
        required=True,
        help="Rooted Newick tree whose leaves match the fixture species IDs.",
    )
    args = parser.parse_args()

    fixture = load_fixture(Path(args.fixture))
    tree = parse_newick(args.tree)
    print(summarize(fixture, tree))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
