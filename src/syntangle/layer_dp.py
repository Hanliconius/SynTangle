from __future__ import annotations

from dataclasses import dataclass
from itertools import permutations
from math import factorial

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
from .orientation_space import legal_orientation_assignments
from .pair_cost import pair_component_crossings


@dataclass(frozen=True)
class ComponentDPDiagnostics:
    component_id: int
    chromosome_count: int
    orientation_assignments: int
    permutation_counts: tuple[tuple[str, int], ...]
    transitions_evaluated: int

    def to_dict(self) -> dict[str, object]:
        return {
            "component_id": self.component_id,
            "chromosome_count": self.chromosome_count,
            "orientation_assignments": self.orientation_assignments,
            "permutation_counts": dict(self.permutation_counts),
            "transitions_evaluated": self.transitions_evaluated,
        }


@dataclass(frozen=True)
class LayerDPResult:
    layout: ExactLayoutResult
    diagnostics: tuple[ComponentDPDiagnostics, ...]

    def to_dict(self) -> dict[str, object]:
        output = self.layout.to_dict()
        output["solver"] = "exact-layer-dynamic-programming"
        output["component_diagnostics"] = [
            item.to_dict() for item in self.diagnostics
        ]
        return output


def _component_data(fixture: Fixture):
    graph = build_incidence_graph(fixture)
    from .component_space import optimization_component_map
    components, component_of = optimization_component_map(fixture)

    return graph, components, component_of


def _best_path_for_orientation(
    fixture: Fixture,
    component_nodes: frozenset[str],
    permutation_options: dict[str, tuple[tuple[ChromosomeRef, ...], ...]],
    orientation: dict[ChromosomeRef, int],
) -> tuple[int, tuple[tuple[ChromosomeRef, ...], ...], int]:
    first_species = fixture.species_ids[0]
    current = {
        order: (0, (order,))
        for order in permutation_options[first_species]
    }
    transitions = 0

    for species1, species2 in zip(
        fixture.species_ids, fixture.species_ids[1:]
    ):
        next_states = {}
        for order2 in permutation_options[species2]:
            best = None
            for order1, (cost_so_far, path_so_far) in current.items():
                pair_cost = pair_component_crossings(
                    fixture,
                    species1,
                    species2,
                    order1,
                    order2,
                    orientation,
                    component_nodes,
                )
                transitions += 1
                path = path_so_far + (order2,)
                path_key = tuple(
                    ref.label for order in path for ref in order
                )
                key = (cost_so_far + pair_cost, path_key)
                if best is None or key < best[0]:
                    best = (key, cost_so_far + pair_cost, path)
            assert best is not None
            next_states[order2] = (best[1], best[2])
        current = next_states

    best_final = None
    for cost, path in current.values():
        path_key = tuple(ref.label for order in path for ref in order)
        key = (cost, path_key)
        if best_final is None or key < best_final[0]:
            best_final = (key, cost, path)

    assert best_final is not None
    return best_final[1], best_final[2], transitions


def exact_optimize_layer_dp(
    fixture: Fixture,
    *,
    orientation_cap_per_component: int = 4096,
    permutation_cap_per_species: int = 40320,
    transition_cap_per_component: int = 5_000_000,
) -> LayerDPResult:
    """Prove the crossing optimum by dynamic programming across species layers.

    For a fixed legal orientation assignment, the crossing objective is a sum
    over adjacent species pairs. Chromosome-order optimization is therefore a
    chain-structured dynamic program rather than a Cartesian product over all
    species permutations.
    """

    initial = initial_layout_state(fixture)
    normalized = canonicalize_component_order(fixture, initial)
    initial_score = score_crossings(fixture, initial)
    normalized_score = score_crossings(fixture, normalized)

    graph, components, component_of = _component_data(fixture)
    chosen_orders: dict[int, dict[str, tuple[ChromosomeRef, ...]]] = {}
    chosen_orientation: dict[int, dict[ChromosomeRef, int]] = {}
    diagnostics: list[ComponentDPDiagnostics] = []
    total_transitions = 0

    for component_id, component_nodes in enumerate(components):
        refs = tuple(
            sorted(
                ref
                for ref in fixture.chromosome_refs
                if component_of[ref] == component_id
            )
        )

        permutation_options = {}
        permutation_counts = []
        for species in fixture.species_ids:
            species_refs = tuple(
                ref for ref in refs if ref.species_id == species
            )
            count = factorial(len(species_refs))
            if count > permutation_cap_per_species:
                raise SearchSpaceTooLarge(
                    f"Component {component_id}, species {species} has "
                    f"{count} chromosome permutations; cap is "
                    f"{permutation_cap_per_species}"
                )
            options = (
                tuple(permutations(species_refs))
                if species_refs
                else ((),)
            )
            permutation_options[species] = options
            permutation_counts.append((species, len(options)))

        orientations = legal_orientation_assignments(
            fixture,
            refs,
            cap=orientation_cap_per_component,
        )

        transition_estimate = len(orientations) * sum(
            len(permutation_options[a]) * len(permutation_options[b])
            for a, b in zip(
                fixture.species_ids, fixture.species_ids[1:]
            )
        )
        if transition_estimate > transition_cap_per_component:
            raise SearchSpaceTooLarge(
                f"Component {component_id} requires up to "
                f"{transition_estimate} DP transitions; cap is "
                f"{transition_cap_per_component}"
            )

        best = None
        component_transitions = 0

        for orientation in orientations:
            cost, path, transitions = _best_path_for_orientation(
                fixture,
                component_nodes,
                permutation_options,
                orientation,
            )
            component_transitions += transitions
            path_key = tuple(
                ref.label for order in path for ref in order
            )
            orientation_key = tuple(
                orientation[ref] for ref in sorted(orientation)
            )
            key = (cost, path_key, orientation_key)

            if best is None or key < best[0]:
                best = (key, path, dict(orientation))

        assert best is not None
        path = best[1]
        chosen_orders[component_id] = {
            species: tuple(order)
            for species, order in zip(fixture.species_ids, path)
            if order
        }
        chosen_orientation[component_id] = best[2]
        total_transitions += component_transitions

        diagnostics.append(
            ComponentDPDiagnostics(
                component_id=component_id,
                chromosome_count=len(refs),
                orientation_assignments=len(orientations),
                permutation_counts=tuple(permutation_counts),
                transitions_evaluated=component_transitions,
            )
        )

    final_order = {}
    for species in fixture.species_ids:
        ordered = []
        for component_id in range(len(components)):
            ordered.extend(
                chosen_orders[component_id].get(species, ())
            )
        final_order[species] = tuple(ordered)

    final_orientation = {}
    for component_id in range(len(components)):
        final_orientation.update(chosen_orientation[component_id])

    optimized = LayoutState(
        chromosome_order=final_order,
        chromosome_orientation=final_orientation,
    )
    optimized_score = score_crossings(fixture, optimized)

    layout = ExactLayoutResult(
        initial_state=initial,
        normalized_state=normalized,
        optimized_state=optimized,
        initial_score=initial_score,
        normalized_score=normalized_score,
        optimized_score=optimized_score,
        states_evaluated=total_transitions,
        optimality_status="proven optimum",
    )

    return LayerDPResult(
        layout=layout,
        diagnostics=tuple(diagnostics),
    )
