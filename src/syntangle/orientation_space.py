from __future__ import annotations

from itertools import product

from .layout import InfeasibleOrientationConstraints, SearchSpaceTooLarge
from .model import ChromosomeRef, Fixture


def legal_orientation_assignments(
    fixture: Fixture,
    refs: tuple[ChromosomeRef, ...],
    *,
    cap: int = 4096,
) -> tuple[dict[ChromosomeRef, int], ...]:
    """Enumerate legal whole-chromosome orientations by GF(2) propagation.

    Connected orientation components each contribute one global reversal bit.
    This avoids filtering all 2^n raw chromosome orientations when hard XOR
    constraints already determine relative states.
    """

    ref_set = set(refs)
    constraints = tuple(
        c
        for c in fixture.orientation_constraints
        if c.a in ref_set or c.b in ref_set
    )
    if any(c.a not in ref_set or c.b not in ref_set for c in constraints):
        raise ValueError(
            "Orientation constraint crosses disconnected incidence components"
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

    n_assignments = 2 ** len(groups)
    if n_assignments > cap:
        raise SearchSpaceTooLarge(
            f"Component has {n_assignments} legal orientation assignments; "
            f"cap is {cap}"
        )

    assignments: list[dict[ChromosomeRef, int]] = []
    for flips in product((0, 1), repeat=len(groups)):
        assignment: dict[ChromosomeRef, int] = {}
        for group, flip in zip(groups, flips):
            for ref in group:
                bit = relative[ref] ^ flip
                assignment[ref] = 1 if bit == 0 else -1
        assignments.append(assignment)

    return tuple(assignments)
