from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from .incidence import build_incidence_graph, chromosome_node_id
from .layout import (
    AmbiguousHomologyError,
    LayoutState,
    _anchor_key,
    _occurrences_by_species_homology,
)
from .model import ChromosomeRef, Fixture


@dataclass(frozen=True)
class OrderDPResult:
    species_id: str
    chromosome_order: tuple[ChromosomeRef, ...]
    variable_crossing_cost: int
    subset_states_evaluated: int

    def to_dict(self) -> dict[str, object]:
        return {
            "species": self.species_id,
            "chromosome_order": [
                ref.chromosome_id for ref in self.chromosome_order
            ],
            "variable_crossing_cost": self.variable_crossing_cost,
            "subset_states_evaluated": self.subset_states_evaluated,
        }


@dataclass(frozen=True)
class PairwiseOrderCosts:
    refs: tuple[ChromosomeRef, ...]
    before_cost: dict[tuple[ChromosomeRef, ChromosomeRef], int]

    def cost(self, before: ChromosomeRef, after: ChromosomeRef) -> int:
        return self.before_cost.get((before, after), 0)


def _neighbor_species(
    fixture: Fixture,
    species_id: str,
) -> tuple[str, ...]:
    index = fixture.species_ids.index(species_id)
    neighbors = []
    if index > 0:
        neighbors.append(fixture.species_ids[index - 1])
    if index + 1 < len(fixture.species_ids):
        neighbors.append(fixture.species_ids[index + 1])
    return tuple(neighbors)


def build_pairwise_order_costs(
    fixture: Fixture,
    state: LayoutState,
    species_id: str,
    component_nodes: frozenset[str],
) -> PairwiseOrderCosts:
    """Build exact pairwise order costs for one species/component.

    Orientations and all neighboring species orders are fixed. The only free
    variables are relative orders among complete chromosomes in this one
    species/component. Crossing contributions involving two different target
    chromosomes are therefore pairwise precedence costs.
    """

    refs = tuple(
        ref
        for ref in state.chromosome_order[species_id]
        if chromosome_node_id(ref) in component_nodes
    )
    if len(refs) < 2:
        return PairwiseOrderCosts(refs=refs, before_cost={})

    occurrence_index = _occurrences_by_species_homology(fixture)
    neighbor_keys: dict[
        str, dict[ChromosomeRef, list[tuple[int, float]]]
    ] = {}

    chromosome_lookup = {
        chromosome.ref: chromosome for chromosome in fixture.chromosomes
    }

    for neighbor in _neighbor_species(fixture, species_id):
        neighbor_rank = {
            ref: rank
            for rank, ref in enumerate(state.chromosome_order[neighbor])
        }
        by_target = {ref: [] for ref in refs}

        shared = sorted(
            set(occurrence_index.get(species_id, {}))
            & set(occurrence_index.get(neighbor, {}))
        )
        for homology_id in shared:
            target_occ = occurrence_index[species_id][homology_id]
            neighbor_occ = occurrence_index[neighbor][homology_id]
            if len(target_occ) != 1 or len(neighbor_occ) != 1:
                raise AmbiguousHomologyError(
                    f"Homology {homology_id!r} is not one-to-one between "
                    f"{species_id} and {neighbor}"
                )

            target_chrom, _ = target_occ[0]
            neighbor_chrom, neighbor_block = neighbor_occ[0]
            if target_chrom.ref not in by_target:
                continue
            if chromosome_node_id(neighbor_chrom.ref) not in component_nodes:
                continue

            key = _anchor_key(
                chromosome_lookup[neighbor_chrom.ref],
                neighbor_block,
                neighbor_rank[neighbor_chrom.ref],
                state.chromosome_orientation[neighbor_chrom.ref],
            )
            by_target[target_chrom.ref].append(key)

        neighbor_keys[neighbor] = by_target

    before_cost: dict[
        tuple[ChromosomeRef, ChromosomeRef], int
    ] = {}

    for a, b in combinations(refs, 2):
        a_before_b = 0
        b_before_a = 0

        for neighbor in neighbor_keys:
            keys_a = neighbor_keys[neighbor][a]
            keys_b = neighbor_keys[neighbor][b]
            for key_a in keys_a:
                for key_b in keys_b:
                    if key_a > key_b:
                        a_before_b += 1
                    elif key_b > key_a:
                        b_before_a += 1

        before_cost[(a, b)] = a_before_b
        before_cost[(b, a)] = b_before_a

    return PairwiseOrderCosts(
        refs=refs,
        before_cost=before_cost,
    )


def solve_order_subset_dp(
    species_id: str,
    costs: PairwiseOrderCosts,
) -> OrderDPResult:
    """Solve one layer's chromosome order exactly in O(n 2^n)."""

    refs = costs.refs
    n = len(refs)
    if n <= 1:
        return OrderDPResult(
            species_id=species_id,
            chromosome_order=refs,
            variable_crossing_cost=0,
            subset_states_evaluated=1,
        )

    best: dict[int, tuple[int, tuple[ChromosomeRef, ...]]] = {
        0: (0, ())
    }
    states_evaluated = 1

    for mask in range(1 << n):
        if mask not in best:
            continue
        base_cost, path = best[mask]

        for index, ref in enumerate(refs):
            bit = 1 << index
            if mask & bit:
                continue

            incremental = sum(
                costs.cost(previous, ref)
                for previous in path
            )
            new_path = path + (ref,)
            new_cost = base_cost + incremental
            new_mask = mask | bit
            states_evaluated += 1

            path_key = tuple(item.label for item in new_path)
            candidate_key = (new_cost, path_key)
            current = best.get(new_mask)

            if current is None:
                best[new_mask] = (new_cost, new_path)
            else:
                current_key = (
                    current[0],
                    tuple(item.label for item in current[1]),
                )
                if candidate_key < current_key:
                    best[new_mask] = (new_cost, new_path)

    final_cost, final_order = best[(1 << n) - 1]
    return OrderDPResult(
        species_id=species_id,
        chromosome_order=final_order,
        variable_crossing_cost=final_cost,
        subset_states_evaluated=states_evaluated,
    )


def optimize_species_component_order(
    fixture: Fixture,
    state: LayoutState,
    species_id: str,
    component_nodes: frozenset[str],
) -> tuple[LayoutState, OrderDPResult]:
    """Return the exact best order for one species/component with neighbors fixed."""

    costs = build_pairwise_order_costs(
        fixture,
        state,
        species_id,
        component_nodes,
    )
    result = solve_order_subset_dp(species_id, costs)
    if len(result.chromosome_order) <= 1:
        return state, result

    component_refs = set(result.chromosome_order)
    positions = [
        index
        for index, ref in enumerate(state.chromosome_order[species_id])
        if ref in component_refs
    ]
    order = dict(state.chromosome_order)
    refs = list(order[species_id])

    for index, ref in zip(positions, result.chromosome_order):
        refs[index] = ref
    order[species_id] = tuple(refs)

    return (
        LayoutState(
            chromosome_order=order,
            chromosome_orientation=dict(state.chromosome_orientation),
        ),
        result,
    )
