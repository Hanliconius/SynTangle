from __future__ import annotations

from dataclasses import dataclass
from itertools import product

from .model import ChromosomeRef, Fixture, OrientationConstraint


@dataclass(frozen=True)
class UnsatisfiedOrientationConstraint:
    index: int
    a: ChromosomeRef
    b: ChromosomeRef
    xor: int

    def to_dict(self) -> dict[str, object]:
        return {
            "index": self.index,
            "a": self.a.label,
            "b": self.b.label,
            "xor": self.xor,
        }


@dataclass(frozen=True)
class OrientationResult:
    balanced: bool
    frustration_index: int | None
    free_bits: int
    assignment: dict[ChromosomeRef, int]
    unsatisfied: tuple[UnsatisfiedOrientationConstraint, ...]
    exact: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "balanced": self.balanced,
            "frustration_index": self.frustration_index,
            "free_bits": self.free_bits,
            "assignment": {ref.label: bit for ref, bit in sorted(self.assignment.items())},
            "unsatisfied": [item.to_dict() for item in self.unsatisfied],
            "exact": self.exact,
        }


def _constraint_adjacency(
    chromosome_refs: tuple[ChromosomeRef, ...],
    constraints: tuple[OrientationConstraint, ...],
) -> dict[ChromosomeRef, list[tuple[ChromosomeRef, int, int]]]:
    adjacency = {ref: [] for ref in chromosome_refs}
    for index, constraint in enumerate(constraints):
        adjacency[constraint.a].append((constraint.b, constraint.xor, index))
        if constraint.a != constraint.b:
            adjacency[constraint.b].append((constraint.a, constraint.xor, index))
    return adjacency


def _connected_components(
    adjacency: dict[ChromosomeRef, list[tuple[ChromosomeRef, int, int]]]
) -> tuple[tuple[ChromosomeRef, ...], ...]:
    unseen = set(adjacency)
    components: list[tuple[ChromosomeRef, ...]] = []
    while unseen:
        root = min(unseen)
        stack = [root]
        component: set[ChromosomeRef] = set()
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            unseen.discard(node)
            for neighbor, _, _ in adjacency[node]:
                if neighbor not in component:
                    stack.append(neighbor)
        components.append(tuple(sorted(component)))
    components.sort(key=lambda c: c[0])
    return tuple(components)


def _propagate(
    chromosome_refs: tuple[ChromosomeRef, ...],
    constraints: tuple[OrientationConstraint, ...],
) -> tuple[dict[ChromosomeRef, int], set[int], int]:
    adjacency = _constraint_adjacency(chromosome_refs, constraints)
    components = _connected_components(adjacency)
    assignment: dict[ChromosomeRef, int] = {}
    contradictions: set[int] = set()

    for component in components:
        root = component[0]
        if root in assignment:
            continue
        assignment[root] = 0
        stack = [root]
        while stack:
            node = stack.pop()
            current = assignment[node]
            for neighbor, xor_value, constraint_index in adjacency[node]:
                expected = current ^ xor_value
                if neighbor not in assignment:
                    assignment[neighbor] = expected
                    stack.append(neighbor)
                elif assignment[neighbor] != expected:
                    contradictions.add(constraint_index)

    return assignment, contradictions, len(components)


def _exact_component_assignment(
    component: tuple[ChromosomeRef, ...],
    constraint_indices: tuple[int, ...],
    constraints: tuple[OrientationConstraint, ...],
) -> tuple[int, dict[ChromosomeRef, int], tuple[int, ...]]:
    # Fix the first chromosome to zero; flipping every bit in a connected
    # component leaves all XOR equations unchanged and removes a redundant
    # factor of two from exhaustive search.
    root = component[0]
    others = component[1:]
    best_cost: int | None = None
    best_assignment: dict[ChromosomeRef, int] | None = None
    best_unsatisfied: tuple[int, ...] | None = None

    for bits in product((0, 1), repeat=len(others)):
        assignment = {root: 0, **dict(zip(others, bits))}
        unsatisfied = tuple(
            index
            for index in constraint_indices
            if (assignment[constraints[index].a] ^ assignment[constraints[index].b])
            != constraints[index].xor
        )
        cost = len(unsatisfied)
        if best_cost is None or cost < best_cost:
            best_cost = cost
            best_assignment = assignment
            best_unsatisfied = unsatisfied
            if cost == 0:
                break

    assert best_cost is not None
    assert best_assignment is not None
    assert best_unsatisfied is not None
    return best_cost, best_assignment, best_unsatisfied


def solve_orientation_constraints(fixture: Fixture, exact_limit: int = 22) -> OrientationResult:
    """Solve whole-chromosome GF(2) orientation equations.

    Balanced components are solved by propagation. If contradictions remain,
    each small connected component is solved exactly for minimum frustration
    after fixing one global-reversal bit. Components larger than exact_limit
    retain a propagated assignment but report frustration_index=None and
    exact=False rather than pretending an optimum was proven.
    """

    refs = fixture.chromosome_refs
    constraints = fixture.orientation_constraints
    propagated, contradictions, free_bits = _propagate(refs, constraints)

    if not contradictions:
        return OrientationResult(
            balanced=True,
            frustration_index=0,
            free_bits=free_bits,
            assignment=propagated,
            unsatisfied=(),
            exact=True,
        )

    adjacency = _constraint_adjacency(refs, constraints)
    components = _connected_components(adjacency)
    exact = True
    total_frustration = 0
    best_assignment: dict[ChromosomeRef, int] = {}
    best_unsatisfied_indices: list[int] = []

    for component in components:
        component_set = set(component)
        constraint_indices = tuple(
            index
            for index, constraint in enumerate(constraints)
            if constraint.a in component_set and constraint.b in component_set
        )

        if not constraint_indices:
            best_assignment[component[0]] = 0
            continue

        if len(component) > exact_limit:
            exact = False
            local_assignment = {ref: propagated[ref] for ref in component}
            best_assignment.update(local_assignment)
            best_unsatisfied_indices.extend(
                index
                for index in constraint_indices
                if (local_assignment[constraints[index].a] ^ local_assignment[constraints[index].b])
                != constraints[index].xor
            )
            continue

        cost, assignment, unsatisfied = _exact_component_assignment(
            component, constraint_indices, constraints
        )
        total_frustration += cost
        best_assignment.update(assignment)
        best_unsatisfied_indices.extend(unsatisfied)

    frustration_index: int | None = total_frustration if exact else None
    unsatisfied = tuple(
        UnsatisfiedOrientationConstraint(
            index=index,
            a=constraints[index].a,
            b=constraints[index].b,
            xor=constraints[index].xor,
        )
        for index in sorted(set(best_unsatisfied_indices))
    )

    return OrientationResult(
        balanced=False,
        frustration_index=frustration_index,
        free_bits=free_bits,
        assignment=best_assignment,
        unsatisfied=unsatisfied,
        exact=exact,
    )
