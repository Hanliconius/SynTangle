from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from math import factorial

from .heuristic import optimize_local_search
from .incidence import build_incidence_graph, chromosome_node_id
from .layout import (
    ExactLayoutResult,
    LayoutState,
    SearchSpaceTooLarge,
    canonicalize_component_order,
    initial_layout_state,
    score_crossings,
)
from .model import ChromosomeRef, Fixture
from .order_dp import build_pairwise_order_costs, solve_order_subset_dp
from .orientation_space import legal_orientation_assignments
from .pair_cost import pair_component_crossings


@dataclass(frozen=True)
class ComponentBranchAndBound:
    component_id: int
    chromosome_count: int
    nodes_evaluated: int
    upper_bound: int
    lower_bound: int
    gap: int
    proven: bool
    orientation_assignments: int

    def to_dict(self) -> dict[str, object]:
        return {
            "component_id": self.component_id,
            "chromosome_count": self.chromosome_count,
            "nodes_evaluated": self.nodes_evaluated,
            "upper_bound": self.upper_bound,
            "lower_bound": self.lower_bound,
            "gap": self.gap,
            "proven": self.proven,
            "orientation_assignments": self.orientation_assignments,
        }


@dataclass(frozen=True)
class BranchAndBoundResult:
    layout: ExactLayoutResult
    lower_bound: int
    upper_bound: int
    gap: int
    components: tuple[ComponentBranchAndBound, ...]

    def to_dict(self) -> dict[str, object]:
        output = self.layout.to_dict()
        output["solver"] = "hard-kernel-branch-and-bound"
        output["lower_bound"] = self.lower_bound
        output["upper_bound"] = self.upper_bound
        output["optimality_gap"] = self.gap
        output["component_diagnostics"] = [
            component.to_dict() for component in self.components
        ]
        return output


def _component_map(
    fixture: Fixture,
) -> tuple[tuple[frozenset[str], ...], dict[ChromosomeRef, int]]:
    graph = build_incidence_graph(fixture)
    components = graph.connected_components()
    component_of: dict[ChromosomeRef, int] = {}
    for component_id, nodes in enumerate(components):
        for ref in fixture.chromosome_refs:
            if chromosome_node_id(ref) in nodes:
                component_of[ref] = component_id
    return components, component_of


def _component_refs_by_species(
    fixture: Fixture,
    component_id: int,
    component_of: dict[ChromosomeRef, int],
) -> dict[str, tuple[ChromosomeRef, ...]]:
    return {
        species: tuple(
            ref
            for ref in fixture.chromosome_refs
            if ref.species_id == species and component_of[ref] == component_id
        )
        for species in fixture.species_ids
    }


def _state_with_component_orders(
    fixture: Fixture,
    base: LayoutState,
    component_id: int,
    component_of: dict[ChromosomeRef, int],
    local_orders: dict[int, tuple[ChromosomeRef, ...]],
    orientation: dict[ChromosomeRef, int],
) -> LayoutState:
    orders = dict(base.chromosome_order)
    for species_index, local_order in local_orders.items():
        species = fixture.species_ids[species_index]
        current = list(orders[species])
        positions = [
            index
            for index, ref in enumerate(current)
            if component_of[ref] == component_id
        ]
        for index, ref in zip(positions, local_order):
            current[index] = ref
        orders[species] = tuple(current)

    full_orientation = dict(base.chromosome_orientation)
    full_orientation.update(orientation)
    return LayoutState(
        chromosome_order=orders,
        chromosome_orientation=full_orientation,
    )


def _component_cost(
    fixture: Fixture,
    state: LayoutState,
    component_nodes: frozenset[str],
) -> int:
    return score_crossings(
        fixture,
        state,
        restrict_component_nodes=component_nodes,
    ).crossings


def _edge_cost(
    fixture: Fixture,
    species_left: str,
    species_right: str,
    order_left: tuple[ChromosomeRef, ...],
    order_right: tuple[ChromosomeRef, ...],
    state: LayoutState,
    component_nodes: frozenset[str],
) -> int:
    return pair_component_crossings(
        fixture,
        species_left,
        species_right,
        order_left,
        order_right,
        state.chromosome_orientation,
        component_nodes,
    )


def _conditional_edge_minimum(
    fixture: Fixture,
    state: LayoutState,
    target_index: int,
    neighbor_index: int,
    component_nodes: frozenset[str],
) -> tuple[int, tuple[ChromosomeRef, ...]]:
    target = fixture.species_ids[target_index]
    neighbor = fixture.species_ids[neighbor_index]
    costs = build_pairwise_order_costs(
        fixture,
        state,
        target,
        component_nodes,
        neighbors=(neighbor,),
    )
    result = solve_order_subset_dp(target, costs)

    orders = dict(state.chromosome_order)
    current = list(orders[target])
    component_refs = set(result.chromosome_order)
    positions = [
        index for index, ref in enumerate(current)
        if ref in component_refs
    ]
    for index, ref in zip(positions, result.chromosome_order):
        current[index] = ref
    orders[target] = tuple(current)
    candidate = LayoutState(
        chromosome_order=orders,
        chromosome_orientation=dict(state.chromosome_orientation),
    )

    if target_index < neighbor_index:
        left_species, right_species = target, neighbor
        left_order = result.chromosome_order
        right_order = tuple(
            ref for ref in candidate.chromosome_order[neighbor]
            if chromosome_node_id(ref) in component_nodes
        )
    else:
        left_species, right_species = neighbor, target
        left_order = tuple(
            ref for ref in candidate.chromosome_order[neighbor]
            if chromosome_node_id(ref) in component_nodes
        )
        right_order = result.chromosome_order

    exact_minimum = _edge_cost(
        fixture,
        left_species,
        right_species,
        left_order,
        right_order,
        candidate,
        component_nodes,
    )
    return exact_minimum, result.chromosome_order


def _frontier_lower_bound(
    fixture: Fixture,
    base: LayoutState,
    component_id: int,
    component_of: dict[ChromosomeRef, int],
    component_nodes: frozenset[str],
    orientation: dict[ChromosomeRef, int],
    local_orders: dict[int, tuple[ChromosomeRef, ...]],
    left: int,
    right: int,
    current_cost: int,
) -> int:
    state = _state_with_component_orders(
        fixture,
        base,
        component_id,
        component_of,
        local_orders,
        orientation,
    )
    bound = current_cost

    if left > 0:
        minimum, _ = _conditional_edge_minimum(
            fixture, state, left - 1, left, component_nodes
        )
        bound += minimum

    if right + 1 < len(fixture.species_ids):
        minimum, _ = _conditional_edge_minimum(
            fixture, state, right + 1, right, component_nodes
        )
        bound += minimum

    return bound


def _ordered_candidate_permutations(
    refs: tuple[ChromosomeRef, ...],
    preferred: tuple[ChromosomeRef, ...] | None,
):
    if len(refs) <= 1:
        yield refs
        return

    seen_preferred = False
    if preferred is not None and set(preferred) == set(refs):
        seen_preferred = True
        yield preferred

    for order in permutations(refs):
        if seen_preferred and order == preferred:
            continue
        yield order


def optimize_branch_and_bound(
    fixture: Fixture,
    *,
    node_cap_per_component: int = 250_000,
    orientation_cap_per_component: int = 4096,
    local_restarts: int = 6,
    seed: int = 1,
) -> BranchAndBoundResult:
    """Exact hard-kernel search with deterministic lower/upper bounds.

    Search begins from the species layer with the fewest chromosomes in each
    incidence component and expands outward. Once one neighboring layer is
    fixed, the exact O(n 2^n) one-layer solver supplies a lower bound for the
    next edge. If the node cap interrupts search, the result is explicitly
    bounded rather than reported as optimal.
    """

    if node_cap_per_component < 1:
        raise ValueError("node_cap_per_component must be at least 1")

    initial = initial_layout_state(fixture)
    normalized = canonicalize_component_order(fixture, initial)
    initial_score = score_crossings(fixture, initial)
    normalized_score = score_crossings(fixture, normalized)

    heuristic = optimize_local_search(
        fixture,
        restarts=local_restarts,
        seed=seed,
    )
    incumbent_state = heuristic.layout.optimized_state

    components, component_of = _component_map(fixture)
    chosen_orders: dict[int, dict[str, tuple[ChromosomeRef, ...]]] = {}
    chosen_orientation: dict[int, dict[ChromosomeRef, int]] = {}
    diagnostics: list[ComponentBranchAndBound] = []

    total_lower = 0
    total_upper = 0
    total_nodes = 0

    for component_id, component_nodes in enumerate(components):
        refs_by_species = _component_refs_by_species(
            fixture, component_id, component_of
        )
        refs = tuple(
            sorted(
                ref for species_refs in refs_by_species.values()
                for ref in species_refs
            )
        )

        incumbent_upper = _component_cost(
            fixture,
            incumbent_state,
            component_nodes,
        )
        best_state = incumbent_state
        nodes = 0
        frontier_lower = float("inf")
        exhausted = True

        if incumbent_upper == 0:
            orientations = ()
        else:
            try:
                orientations = legal_orientation_assignments(
                    fixture,
                    refs,
                    cap=orientation_cap_per_component,
                )
            except SearchSpaceTooLarge:
                orientations = ()
                exhausted = False
                frontier_lower = 0

        if incumbent_upper == 0:
            exhausted = True
            frontier_lower = 0
        elif orientations:
            nonempty_indices = [
                index
                for index, species in enumerate(fixture.species_ids)
                if refs_by_species[species]
            ]
            pivot = min(
                nonempty_indices,
                key=lambda index: (
                    len(refs_by_species[fixture.species_ids[index]]),
                    abs(index - (len(fixture.species_ids) - 1) / 2),
                    index,
                ),
            )

            stop_all = False

            for orientation_index, orientation in enumerate(orientations):
                if stop_all:
                    exhausted = False
                    frontier_lower = min(frontier_lower, 0)
                    break

                pivot_refs = refs_by_species[
                    fixture.species_ids[pivot]
                ]
                preferred_pivot = tuple(
                    ref
                    for ref in incumbent_state.chromosome_order[
                        fixture.species_ids[pivot]
                    ]
                    if component_of[ref] == component_id
                )

                for pivot_order in _ordered_candidate_permutations(
                    pivot_refs, preferred_pivot
                ):
                    if nodes >= node_cap_per_component:
                        exhausted = False
                        frontier_lower = min(frontier_lower, 0)
                        stop_all = True
                        break

                    local_orders = {pivot: pivot_order}
                    root_state = _state_with_component_orders(
                        fixture,
                        normalized,
                        component_id,
                        component_of,
                        local_orders,
                        orientation,
                    )
                    root_bound = _frontier_lower_bound(
                        fixture,
                        normalized,
                        component_id,
                        component_of,
                        component_nodes,
                        orientation,
                        local_orders,
                        pivot,
                        pivot,
                        0,
                    )
                    nodes += 1
                    if root_bound >= incumbent_upper:
                        continue

                    def search(
                        left: int,
                        right: int,
                        orders: dict[int, tuple[ChromosomeRef, ...]],
                        current_cost: int,
                        node_bound: int,
                    ) -> None:
                        nonlocal nodes, incumbent_upper, best_state
                        nonlocal exhausted, frontier_lower, stop_all

                        if stop_all:
                            frontier_lower = min(frontier_lower, node_bound)
                            return
                        if current_cost >= incumbent_upper:
                            return
                        if left == 0 and right == len(fixture.species_ids) - 1:
                            candidate = _state_with_component_orders(
                                fixture,
                                normalized,
                                component_id,
                                component_of,
                                orders,
                                orientation,
                            )
                            exact_cost = _component_cost(
                                fixture, candidate, component_nodes
                            )
                            if exact_cost < incumbent_upper:
                                incumbent_upper = exact_cost
                                best_state = candidate
                            elif exact_cost == incumbent_upper:
                                candidate_key = tuple(
                                    ref.label
                                    for species in fixture.species_ids
                                    for ref in candidate.chromosome_order[species]
                                    if component_of[ref] == component_id
                                )
                                best_key = tuple(
                                    ref.label
                                    for species in fixture.species_ids
                                    for ref in best_state.chromosome_order[species]
                                    if component_of[ref] == component_id
                                )
                                if candidate_key < best_key:
                                    best_state = candidate
                            return

                        candidates = []
                        if left > 0:
                            candidates.append(left - 1)
                        if right + 1 < len(fixture.species_ids):
                            candidates.append(right + 1)

                        target = min(
                            candidates,
                            key=lambda index: (
                                factorial(
                                    len(
                                        refs_by_species[
                                            fixture.species_ids[index]
                                        ]
                                    )
                                ),
                                index,
                            ),
                        )
                        neighbor = left if target < left else right

                        current_state = _state_with_component_orders(
                            fixture,
                            normalized,
                            component_id,
                            component_of,
                            orders,
                            orientation,
                        )
                        _, preferred = _conditional_edge_minimum(
                            fixture,
                            current_state,
                            target,
                            neighbor,
                            component_nodes,
                        )
                        target_refs = refs_by_species[
                            fixture.species_ids[target]
                        ]

                        for target_order in _ordered_candidate_permutations(
                            target_refs, preferred
                        ):
                            if nodes >= node_cap_per_component:
                                exhausted = False
                                frontier_lower = min(
                                    frontier_lower, node_bound
                                )
                                stop_all = True
                                return

                            next_orders = dict(orders)
                            next_orders[target] = target_order
                            next_state = _state_with_component_orders(
                                fixture,
                                normalized,
                                component_id,
                                component_of,
                                next_orders,
                                orientation,
                            )

                            if target < neighbor:
                                edge = _edge_cost(
                                    fixture,
                                    fixture.species_ids[target],
                                    fixture.species_ids[neighbor],
                                    target_order,
                                    orders[neighbor],
                                    next_state,
                                    component_nodes,
                                )
                                next_left, next_right = target, right
                            else:
                                edge = _edge_cost(
                                    fixture,
                                    fixture.species_ids[neighbor],
                                    fixture.species_ids[target],
                                    orders[neighbor],
                                    target_order,
                                    next_state,
                                    component_nodes,
                                )
                                next_left, next_right = left, target

                            next_cost = current_cost + edge
                            if next_cost >= incumbent_upper:
                                continue

                            next_bound = _frontier_lower_bound(
                                fixture,
                                normalized,
                                component_id,
                                component_of,
                                component_nodes,
                                orientation,
                                next_orders,
                                next_left,
                                next_right,
                                next_cost,
                            )
                            nodes += 1

                            if next_bound >= incumbent_upper:
                                continue

                            search(
                                next_left,
                                next_right,
                                next_orders,
                                next_cost,
                                next_bound,
                            )
                            if stop_all:
                                return

                    search(
                        pivot,
                        pivot,
                        local_orders,
                        0,
                        root_bound,
                    )
                    if stop_all:
                        break

        component_lower = (
            incumbent_upper
            if exhausted
            else int(min(frontier_lower, incumbent_upper))
        )
        proven = exhausted or component_lower == incumbent_upper

        local_best_orders: dict[str, tuple[ChromosomeRef, ...]] = {}
        for species in fixture.species_ids:
            local_best_orders[species] = tuple(
                ref
                for ref in best_state.chromosome_order[species]
                if component_of[ref] == component_id
            )

        chosen_orders[component_id] = local_best_orders
        chosen_orientation[component_id] = {
            ref: best_state.chromosome_orientation[ref]
            for ref in refs
        }

        diagnostics.append(
            ComponentBranchAndBound(
                component_id=component_id,
                chromosome_count=len(refs),
                nodes_evaluated=nodes,
                upper_bound=incumbent_upper,
                lower_bound=component_lower,
                gap=incumbent_upper - component_lower,
                proven=proven,
                orientation_assignments=len(orientations),
            )
        )
        total_lower += component_lower
        total_upper += incumbent_upper
        total_nodes += nodes

    final_order: dict[str, tuple[ChromosomeRef, ...]] = {}
    for species in fixture.species_ids:
        ordered: list[ChromosomeRef] = []
        for component_id in range(len(components)):
            ordered.extend(chosen_orders[component_id][species])
        final_order[species] = tuple(ordered)

    final_orientation: dict[ChromosomeRef, int] = {}
    for component_id in range(len(components)):
        final_orientation.update(chosen_orientation[component_id])

    optimized = LayoutState(
        chromosome_order=final_order,
        chromosome_orientation=final_orientation,
    )
    optimized_score = score_crossings(fixture, optimized)
    total_upper = optimized_score.crossings
    total_lower = min(total_lower, total_upper)
    gap = total_upper - total_lower
    status = "proven optimum" if gap == 0 else "bounded best known"

    layout = ExactLayoutResult(
        initial_state=initial,
        normalized_state=normalized,
        optimized_state=optimized,
        initial_score=initial_score,
        normalized_score=normalized_score,
        optimized_score=optimized_score,
        states_evaluated=total_nodes,
        optimality_status=status,
    )

    return BranchAndBoundResult(
        layout=layout,
        lower_bound=total_lower,
        upper_bound=total_upper,
        gap=gap,
        components=tuple(diagnostics),
    )
