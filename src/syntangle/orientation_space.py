from __future__ import annotations

from dataclasses import dataclass
from itertools import product

from .layout import InfeasibleOrientationConstraints, SearchSpaceTooLarge
from .model import ChromosomeRef, Fixture


@dataclass(frozen=True)
class OrientationBasis:
    base_assignment: dict[ChromosomeRef, int]
    free_flip_groups: tuple[tuple[ChromosomeRef, ...], ...]

    @property
    def assignment_count(self) -> int:
        return 2 ** len(self.free_flip_groups)


def orientation_basis(
    fixture: Fixture,
    refs: tuple[ChromosomeRef, ...],
) -> OrientationBasis:
    """Compactly represent every orientation satisfying hard GF(2) constraints."""

    ref_set = set(refs)
    constraints = tuple(
        constraint
        for constraint in fixture.orientation_constraints
        if constraint.a in ref_set or constraint.b in ref_set
    )
    if any(
        constraint.a not in ref_set or constraint.b not in ref_set
        for constraint in constraints
    ):
        raise ValueError(
            "Orientation constraint crosses the supplied chromosome subset"
        )

    adjacency: dict[ChromosomeRef, list[tuple[ChromosomeRef, int]]] = {
        ref: [] for ref in refs
    }
    for constraint in constraints:
        adjacency[constraint.a].append((constraint.b, constraint.xor))
        if constraint.a != constraint.b:
            adjacency[constraint.b].append((constraint.a, constraint.xor))

    relative: dict[ChromosomeRef, int] = {}
    groups: list[tuple[ChromosomeRef, ...]] = []

    for root in sorted(refs):
        if root in relative:
            continue
        relative[root] = 0
        stack = [root]
        group: set[ChromosomeRef] = set()

        while stack:
            node = stack.pop()
            if node in group:
                continue
            group.add(node)
            for neighbor, xor_value in adjacency[node]:
                expected = relative[node] ^ xor_value
                if neighbor not in relative:
                    relative[neighbor] = expected
                    stack.append(neighbor)
                elif relative[neighbor] != expected:
                    raise InfeasibleOrientationConstraints(
                        "No whole-chromosome orientation satisfies all hard "
                        "orientation constraints"
                    )

        groups.append(tuple(sorted(group)))

    base_assignment = {
        ref: (1 if relative[ref] == 0 else -1)
        for ref in refs
    }
    return OrientationBasis(
        base_assignment=base_assignment,
        free_flip_groups=tuple(groups),
    )


def legal_orientation_assignments(
    fixture: Fixture,
    refs: tuple[ChromosomeRef, ...],
    *,
    cap: int = 4096,
) -> tuple[dict[ChromosomeRef, int], ...]:
    """Enumerate legal orientations from the compact propagated basis."""

    basis = orientation_basis(fixture, refs)
    if basis.assignment_count > cap:
        raise SearchSpaceTooLarge(
            f"Component has {basis.assignment_count} legal orientation "
            f"assignments; cap is {cap}"
        )

    assignments: list[dict[ChromosomeRef, int]] = []
    for flips in product((0, 1), repeat=len(basis.free_flip_groups)):
        assignment = dict(basis.base_assignment)
        for group, flip in zip(basis.free_flip_groups, flips):
            if not flip:
                continue
            for ref in group:
                assignment[ref] *= -1
        assignments.append(assignment)

    return tuple(assignments)
