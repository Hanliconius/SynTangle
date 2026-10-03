from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import heapq
import multiprocessing as mp
from itertools import permutations
from math import factorial

from .bounds import RelaxedCrossingBound, build_relaxed_crossing_bound
from .heuristic import optimize_local_search
from .incidence import build_incidence_graph, chromosome_node_id
from .layout import (
    ExactLayoutResult,
    LayoutState,
    canonicalize_component_order,
    initial_layout_state,
    score_crossings,
)
from .model import ChromosomeRef, Fixture
from .order_dp import build_pairwise_order_costs, solve_order_subset_dp
from .orientation_space import OrientationBasis, orientation_basis
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
    free_orientation_groups: int = 0
    implicit_orientation_states: int = 0
    root_lower_bound: int = 0
    root_unresolved_orientation_groups: int = 0
    orientation_nodes_evaluated: int = 0
    orientation_branches_pruned: int = 0
    orientation_groups_forced: int = 0
    order_nodes_evaluated: int = 0
    memo_hits: int = 0

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
            "free_orientation_groups": self.free_orientation_groups,
            "implicit_orientation_states": self.implicit_orientation_states,
            "root_lower_bound": self.root_lower_bound,
            "root_unresolved_orientation_groups": (
                self.root_unresolved_orientation_groups
            ),
            "orientation_nodes_evaluated": self.orientation_nodes_evaluated,
            "orientation_branches_pruned": self.orientation_branches_pruned,
            "orientation_groups_forced": self.orientation_groups_forced,
            "order_nodes_evaluated": self.order_nodes_evaluated,
            "memo_hits": self.memo_hits,
        }


@dataclass(frozen=True)
class BranchAndBoundResult:
    layout: ExactLayoutResult
    lower_bound: int
    upper_bound: int
    gap: int
    components: tuple[ComponentBranchAndBound, ...]
    component_workers: int = 1

    def to_dict(self) -> dict[str, object]:
        output = self.layout.to_dict()
        output["solver"] = "monotone-component-branch-and-bound"
        output["lower_bound"] = self.lower_bound
        output["upper_bound"] = self.upper_bound
        output["optimality_gap"] = self.gap
        output["component_diagnostics"] = [
            component.to_dict() for component in self.components
        ]
        output["component_workers"] = self.component_workers
        return output


@dataclass
class _Budget:
    cap: int
    used: int = 0

    def consume(self) -> bool:
        if self.used >= self.cap:
            return False
        self.used += 1
        return True

    @property
    def exhausted(self) -> bool:
        return self.used >= self.cap


@dataclass(frozen=True)
class _OrderSearchResult:
    best_state: LayoutState
    upper_bound: int
    lower_bound: int
    proven: bool
    nodes: int
    memo_hits: int


@dataclass(frozen=True)
class _ReducedOrientation:
    bits: tuple[int | None, ...]
    lower_bound: int
    forced: int
    pruned: bool


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
            if ref.species_id == species
            and component_of[ref] == component_id
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
    *,
    cache: dict,
) -> tuple[int, tuple[ChromosomeRef, ...], bool]:
    target = fixture.species_ids[target_index]
    neighbor = fixture.species_ids[neighbor_index]
    neighbor_order = tuple(
        ref
        for ref in state.chromosome_order[neighbor]
        if chromosome_node_id(ref) in component_nodes
    )
    orientation_key = tuple(
        (ref.label, state.chromosome_orientation[ref])
        for ref in sorted(state.chromosome_orientation)
        if chromosome_node_id(ref) in component_nodes
    )
    key = (
        orientation_key,
        target_index,
        neighbor_index,
        tuple(ref.label for ref in neighbor_order),
    )

    cached = cache.get(key)
    if cached is not None:
        return cached[0], cached[1], True

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
        index
        for index, ref in enumerate(current)
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
        right_order = neighbor_order
    else:
        left_species, right_species = neighbor, target
        left_order = neighbor_order
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
    cache[key] = (exact_minimum, result.chromosome_order)
    return exact_minimum, result.chromosome_order, False


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
    relaxed_lower_bound: int,
    conditional_cache: dict,
) -> tuple[int, int]:
    state = _state_with_component_orders(
        fixture,
        base,
        component_id,
        component_of,
        local_orders,
        orientation,
    )
    bound = current_cost
    memo_hits = 0

    if left > 0:
        minimum, _, hit = _conditional_edge_minimum(
            fixture,
            state,
            left - 1,
            left,
            component_nodes,
            cache=conditional_cache,
        )
        bound += minimum
        memo_hits += int(hit)

    if right + 1 < len(fixture.species_ids):
        minimum, _, hit = _conditional_edge_minimum(
            fixture,
            state,
            right + 1,
            right,
            component_nodes,
            cache=conditional_cache,
        )
        bound += minimum
        memo_hits += int(hit)

    return max(bound, relaxed_lower_bound), memo_hits


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


def _complete_orientation(
    basis: OrientationBasis,
    bits: tuple[int | None, ...],
) -> dict[ChromosomeRef, int]:
    if any(bit is None for bit in bits):
        raise ValueError("Cannot materialize a partial orientation assignment")

    orientation = dict(basis.base_assignment)
    for group, bit in zip(basis.free_flip_groups, bits):
        if bit:
            for ref in group:
                orientation[ref] *= -1
    return orientation


def _bit_key(bits: tuple[int | None, ...]) -> tuple[int, ...]:
    return tuple(-1 if bit is None else bit for bit in bits)


def _reduce_orientation(
    bits: tuple[int | None, ...],
    *,
    incumbent_upper: int,
    bound: RelaxedCrossingBound,
    cache: dict[tuple[int, ...], int],
) -> tuple[_ReducedOrientation, int]:
    working = list(bits)
    forced = 0
    memo_hits = 0

    def lower(candidate: tuple[int | None, ...]) -> int:
        nonlocal memo_hits
        key = _bit_key(candidate)
        if key in cache:
            memo_hits += 1
            return cache[key]
        value = bound.lower_bound(candidate)
        cache[key] = value
        return value

    while True:
        current = tuple(working)
        current_bound = lower(current)
        if current_bound >= incumbent_upper:
            return (
                _ReducedOrientation(
                    bits=current,
                    lower_bound=current_bound,
                    forced=forced,
                    pruned=True,
                ),
                memo_hits,
            )

        changed = False
        for index, value in enumerate(working):
            if value is not None:
                continue

            alternatives: list[int] = []
            for bit in (0, 1):
                trial = list(working)
                trial[index] = bit
                alternatives.append(lower(tuple(trial)))

            if (
                alternatives[0] >= incumbent_upper
                and alternatives[1] >= incumbent_upper
            ):
                return (
                    _ReducedOrientation(
                        bits=tuple(working),
                        lower_bound=min(alternatives),
                        forced=forced,
                        pruned=True,
                    ),
                    memo_hits,
                )

            if alternatives[0] >= incumbent_upper:
                working[index] = 1
                forced += 1
                changed = True
                break

            if alternatives[1] >= incumbent_upper:
                working[index] = 0
                forced += 1
                changed = True
                break

        if not changed:
            final = tuple(working)
            return (
                _ReducedOrientation(
                    bits=final,
                    lower_bound=lower(final),
                    forced=forced,
                    pruned=False,
                ),
                memo_hits,
            )


def _choose_orientation_branch(
    bits: tuple[int | None, ...],
    bound: RelaxedCrossingBound,
    cache: dict[tuple[int, ...], int],
) -> tuple[int, int]:
    best_index = -1
    best_key: tuple[int, int, int, int] | None = None
    memo_hits = 0

    for index, value in enumerate(bits):
        if value is not None:
            continue

        child_bounds: list[int] = []
        for bit in (0, 1):
            trial = list(bits)
            trial[index] = bit
            key = _bit_key(tuple(trial))
            if key in cache:
                memo_hits += 1
                value_bound = cache[key]
            else:
                value_bound = bound.lower_bound(tuple(trial))
                cache[key] = value_bound
            child_bounds.append(value_bound)

        key = (
            min(child_bounds),
            max(child_bounds),
            abs(child_bounds[0] - child_bounds[1]),
            -index,
        )
        if best_key is None or key > best_key:
            best_key = key
            best_index = index

    if best_index < 0:
        raise ValueError("No unresolved orientation group remains")
    return best_index, memo_hits


def _search_orders_for_orientation(
    fixture: Fixture,
    normalized: LayoutState,
    component_id: int,
    component_of: dict[ChromosomeRef, int],
    component_nodes: frozenset[str],
    refs_by_species: dict[str, tuple[ChromosomeRef, ...]],
    orientation: dict[ChromosomeRef, int],
    relaxed_lower_bound: int,
    incumbent_state: LayoutState,
    incumbent_upper: int,
    budget: _Budget,
    conditional_cache: dict,
) -> _OrderSearchResult:
    if relaxed_lower_bound >= incumbent_upper:
        return _OrderSearchResult(
            best_state=incumbent_state,
            upper_bound=incumbent_upper,
            lower_bound=relaxed_lower_bound,
            proven=True,
            nodes=0,
            memo_hits=0,
        )

    nonempty_indices = [
        index
        for index, species in enumerate(fixture.species_ids)
        if refs_by_species[species]
    ]
    if not nonempty_indices:
        return _OrderSearchResult(
            best_state=incumbent_state,
            upper_bound=incumbent_upper,
            lower_bound=0,
            proven=True,
            nodes=0,
            memo_hits=0,
        )

    pivot = min(
        nonempty_indices,
        key=lambda index: (
            len(refs_by_species[fixture.species_ids[index]]),
            abs(index - (len(fixture.species_ids) - 1) / 2),
            index,
        ),
    )

    best_state = incumbent_state
    best_upper = incumbent_upper
    nodes = 0
    memo_hits = 0
    exhausted = True
    open_lower = float("inf")
    stop_all = False

    pivot_refs = refs_by_species[fixture.species_ids[pivot]]
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
        local_orders = {pivot: pivot_order}
        root_bound, hits = _frontier_lower_bound(
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
            relaxed_lower_bound,
            conditional_cache,
        )
        memo_hits += hits

        if root_bound >= best_upper:
            continue

        if not budget.consume():
            exhausted = False
            open_lower = min(open_lower, relaxed_lower_bound)
            break
        nodes += 1

        def search(
            left: int,
            right: int,
            orders: dict[int, tuple[ChromosomeRef, ...]],
            current_cost: int,
            node_bound: int,
        ) -> None:
            nonlocal nodes, memo_hits, best_upper, best_state
            nonlocal exhausted, open_lower, stop_all

            if stop_all:
                open_lower = min(open_lower, node_bound)
                return
            if node_bound >= best_upper or current_cost >= best_upper:
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
                if exact_cost < best_upper:
                    best_upper = exact_cost
                    best_state = candidate
                elif exact_cost == best_upper:
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
            _, preferred, hit = _conditional_edge_minimum(
                fixture,
                current_state,
                target,
                neighbor,
                component_nodes,
                cache=conditional_cache,
            )
            memo_hits += int(hit)

            target_refs = refs_by_species[
                fixture.species_ids[target]
            ]

            for target_order in _ordered_candidate_permutations(
                target_refs, preferred
            ):
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
                if next_cost >= best_upper:
                    continue

                next_bound, hits = _frontier_lower_bound(
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
                    relaxed_lower_bound,
                    conditional_cache,
                )
                memo_hits += hits
                if next_bound >= best_upper:
                    continue

                if not budget.consume():
                    exhausted = False
                    open_lower = min(open_lower, node_bound)
                    stop_all = True
                    return
                nodes += 1

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

    if exhausted:
        lower = best_upper
    else:
        lower = int(min(open_lower, relaxed_lower_bound, best_upper))

    return _OrderSearchResult(
        best_state=best_state,
        upper_bound=best_upper,
        lower_bound=lower,
        proven=exhausted or lower == best_upper,
        nodes=nodes,
        memo_hits=memo_hits,
    )


@dataclass(frozen=True)
class _ComponentSolveResult:
    component_id: int
    orders: dict[str, tuple[ChromosomeRef, ...]]
    orientation: dict[ChromosomeRef, int]
    lower_bound: int
    nodes_evaluated: int
    diagnostic: ComponentBranchAndBound


def _solve_branch_component(
    fixture: Fixture,
    normalized: LayoutState,
    incumbent_state: LayoutState,
    component_id: int,
    component_nodes: frozenset[str],
    component_of: dict[ChromosomeRef, int],
    node_cap_per_component: int,
) -> _ComponentSolveResult:
    """Solve one exact incidence component.

    This function is deliberately self-contained and pickle-safe so independent
    incidence components can be dispatched to separate worker processes.
    """

    refs_by_species = _component_refs_by_species(
        fixture, component_id, component_of
    )
    refs = tuple(
        sorted(
            ref
            for species_refs in refs_by_species.values()
            for ref in species_refs
        )
    )

    basis = orientation_basis(fixture, refs)
    relaxed_bound = build_relaxed_crossing_bound(
        fixture, component_nodes, basis
    )
    root_bits: tuple[int | None, ...] = tuple(
        None for _ in basis.free_flip_groups
    )
    root_lower_bound = relaxed_bound.lower_bound(root_bits)

    incumbent_upper = _component_cost(
        fixture,
        incumbent_state,
        component_nodes,
    )
    if root_lower_bound > incumbent_upper:
        raise AssertionError(
            "Relaxed lower bound exceeds a feasible incumbent; "
            "the branch-and-bound bound is invalid"
        )

    best_state = incumbent_state
    budget = _Budget(node_cap_per_component)
    lower_cache: dict[tuple[int, ...], int] = {
        _bit_key(root_bits): root_lower_bound
    }
    conditional_cache: dict = {}

    orientation_nodes = 0
    orientation_leaves = 0
    orientation_pruned = 0
    orientation_forced = 0
    order_nodes = 0
    memo_hits = 0
    unresolved_order_lower: int | None = None

    if incumbent_upper == 0 or root_lower_bound == incumbent_upper:
        component_lower = incumbent_upper
        proven = True
        root_reduced = _ReducedOrientation(
            bits=root_bits,
            lower_bound=root_lower_bound,
            forced=0,
            pruned=True,
        )
    else:
        root_reduced, hits = _reduce_orientation(
            root_bits,
            incumbent_upper=incumbent_upper,
            bound=relaxed_bound,
            cache=lower_cache,
        )
        memo_hits += hits
        orientation_forced += root_reduced.forced

        if root_reduced.pruned:
            component_lower = incumbent_upper
            proven = True
        else:
            root_unresolved = sum(
                bit is None for bit in root_reduced.bits
            )
            heap: list[
                tuple[
                    int,
                    int,
                    tuple[int, ...],
                    tuple[int | None, ...],
                ]
            ] = []
            heapq.heappush(
                heap,
                (
                    root_reduced.lower_bound,
                    root_unresolved,
                    _bit_key(root_reduced.bits),
                    root_reduced.bits,
                ),
            )
            seen = {_bit_key(root_reduced.bits)}

            while heap:
                if heap[0][0] >= incumbent_upper:
                    heap.clear()
                    break
                if budget.exhausted:
                    break

                node_lower, _, _, bits = heapq.heappop(heap)

                refreshed, hits = _reduce_orientation(
                    bits,
                    incumbent_upper=incumbent_upper,
                    bound=relaxed_bound,
                    cache=lower_cache,
                )
                memo_hits += hits
                orientation_forced += refreshed.forced
                if refreshed.pruned:
                    orientation_pruned += 1
                    continue
                bits = refreshed.bits
                node_lower = refreshed.lower_bound

                if node_lower >= incumbent_upper:
                    orientation_pruned += 1
                    continue

                if not budget.consume():
                    heapq.heappush(
                        heap,
                        (
                            node_lower,
                            sum(bit is None for bit in bits),
                            _bit_key(bits),
                            bits,
                        ),
                    )
                    break

                orientation_nodes += 1

                if all(bit is not None for bit in bits):
                    orientation_leaves += 1
                    orientation = _complete_orientation(
                        basis, bits
                    )
                    fixed_lower = relaxed_bound.lower_bound(bits)
                    order_result = _search_orders_for_orientation(
                        fixture,
                        normalized,
                        component_id,
                        component_of,
                        component_nodes,
                        refs_by_species,
                        orientation,
                        fixed_lower,
                        best_state,
                        incumbent_upper,
                        budget,
                        conditional_cache,
                    )
                    order_nodes += order_result.nodes
                    memo_hits += order_result.memo_hits

                    if order_result.upper_bound < incumbent_upper:
                        incumbent_upper = order_result.upper_bound
                        best_state = order_result.best_state
                    elif (
                        order_result.upper_bound == incumbent_upper
                        and order_result.best_state is not best_state
                    ):
                        best_state = order_result.best_state

                    if not order_result.proven:
                        if unresolved_order_lower is None:
                            unresolved_order_lower = (
                                order_result.lower_bound
                            )
                        else:
                            unresolved_order_lower = min(
                                unresolved_order_lower,
                                order_result.lower_bound,
                            )
                        break
                    continue

                branch_index, hits = _choose_orientation_branch(
                    bits,
                    relaxed_bound,
                    lower_cache,
                )
                memo_hits += hits

                for bit in (0, 1):
                    child = list(bits)
                    child[branch_index] = bit
                    reduced, hits = _reduce_orientation(
                        tuple(child),
                        incumbent_upper=incumbent_upper,
                        bound=relaxed_bound,
                        cache=lower_cache,
                    )
                    memo_hits += hits
                    orientation_forced += reduced.forced

                    if (
                        reduced.pruned
                        or reduced.lower_bound >= incumbent_upper
                    ):
                        orientation_pruned += 1
                        continue

                    key = _bit_key(reduced.bits)
                    if key in seen:
                        memo_hits += 1
                        continue
                    seen.add(key)
                    heapq.heappush(
                        heap,
                        (
                            reduced.lower_bound,
                            sum(
                                value is None
                                for value in reduced.bits
                            ),
                            key,
                            reduced.bits,
                        ),
                    )

            open_bounds: list[int] = []
            if heap:
                open_bounds.append(heap[0][0])
            if unresolved_order_lower is not None:
                open_bounds.append(unresolved_order_lower)

            if open_bounds:
                component_lower = min(
                    min(open_bounds),
                    incumbent_upper,
                )
                proven = component_lower == incumbent_upper
            else:
                component_lower = incumbent_upper
                proven = True

    root_unresolved_groups = sum(
        bit is None for bit in root_reduced.bits
    )

    local_best_orders: dict[str, tuple[ChromosomeRef, ...]] = {}
    for species in fixture.species_ids:
        local_best_orders[species] = tuple(
            ref
            for ref in best_state.chromosome_order[species]
            if component_of[ref] == component_id
        )

    local_orientation = {
        ref: best_state.chromosome_orientation[ref]
        for ref in refs
    }
    component_nodes_evaluated = orientation_nodes + order_nodes
    diagnostic = ComponentBranchAndBound(
        component_id=component_id,
        chromosome_count=len(refs),
        nodes_evaluated=component_nodes_evaluated,
        upper_bound=incumbent_upper,
        lower_bound=component_lower,
        gap=incumbent_upper - component_lower,
        proven=proven,
        orientation_assignments=orientation_leaves,
        free_orientation_groups=len(basis.free_flip_groups),
        implicit_orientation_states=basis.assignment_count,
        root_lower_bound=root_lower_bound,
        root_unresolved_orientation_groups=root_unresolved_groups,
        orientation_nodes_evaluated=orientation_nodes,
        orientation_branches_pruned=orientation_pruned,
        orientation_groups_forced=orientation_forced,
        order_nodes_evaluated=order_nodes,
        memo_hits=memo_hits,
    )
    return _ComponentSolveResult(
        component_id=component_id,
        orders=local_best_orders,
        orientation=local_orientation,
        lower_bound=component_lower,
        nodes_evaluated=component_nodes_evaluated,
        diagnostic=diagnostic,
    )


def optimize_branch_and_bound(
    fixture: Fixture,
    *,
    node_cap_per_component: int = 250_000,
    orientation_cap_per_component: int = 4096,
    local_restarts: int = 6,
    seed: int = 1,
    component_workers: int = 1,
) -> BranchAndBoundResult:
    """Monotone exact/bounded search over residual orientation and order choices.

    Independent chromosome-homology incidence components are exact additive
    factors of the current crossing objective. When component_workers > 1,
    those components are solved in separate processes and then recombined
    deterministically. Within each component, hard orientation equations remain
    in their compact GF(2) basis and Stage-15 monotone pruning is unchanged.

    The orientation_cap_per_component argument is retained for API
    compatibility with earlier releases. It no longer caps branch-and-bound
    orientation states; the shared node cap controls actual search work.
    """

    if node_cap_per_component < 1:
        raise ValueError("node_cap_per_component must be at least 1")
    if orientation_cap_per_component < 1:
        raise ValueError("orientation_cap_per_component must be at least 1")
    if component_workers < 1:
        raise ValueError("component_workers must be at least 1")

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
    worker_count = min(component_workers, max(1, len(components)))

    if worker_count > 1 and len(components) > 1:
        with ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=mp.get_context("spawn"),
        ) as executor:
            futures = [
                executor.submit(
                    _solve_branch_component,
                    fixture,
                    normalized,
                    incumbent_state,
                    component_id,
                    component_nodes,
                    component_of,
                    node_cap_per_component,
                )
                for component_id, component_nodes in enumerate(components)
            ]
            component_results = [future.result() for future in futures]
    else:
        component_results = [
            _solve_branch_component(
                fixture,
                normalized,
                incumbent_state,
                component_id,
                component_nodes,
                component_of,
                node_cap_per_component,
            )
            for component_id, component_nodes in enumerate(components)
        ]

    component_results.sort(key=lambda item: item.component_id)
    chosen_orders = {
        item.component_id: item.orders for item in component_results
    }
    chosen_orientation = {
        item.component_id: item.orientation for item in component_results
    }
    diagnostics = tuple(
        item.diagnostic for item in component_results
    )
    total_lower = sum(
        item.lower_bound for item in component_results
    )
    total_nodes = sum(
        item.nodes_evaluated for item in component_results
    )

    final_order: dict[str, tuple[ChromosomeRef, ...]] = {}
    for species in fixture.species_ids:
        ordered: list[ChromosomeRef] = []
        for component_id in range(len(components)):
            ordered.extend(chosen_orders[component_id][species])
        final_order[species] = tuple(ordered)

    final_orientation: dict[ChromosomeRef, int] = {}
    for component_id in range(len(components)):
        final_orientation.update(
            chosen_orientation[component_id]
        )

    optimized = LayoutState(
        chromosome_order=final_order,
        chromosome_orientation=final_orientation,
    )
    optimized_score = score_crossings(fixture, optimized)
    total_upper = optimized_score.crossings
    total_lower = min(total_lower, total_upper)
    gap = total_upper - total_lower
    status = (
        "proven optimum"
        if gap == 0
        else "bounded best known"
    )

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
        components=diagnostics,
        component_workers=worker_count,
    )
