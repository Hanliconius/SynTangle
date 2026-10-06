from __future__ import annotations

from dataclasses import dataclass
import random
from typing import Callable

from .crossing_cache import CrossingCostCache

from .incidence import build_incidence_graph, chromosome_node_id
from .layout import (
    ExactLayoutResult,
    LayoutState,
    canonicalize_component_order,
    initial_layout_state,
    score_crossings,
)
from .model import ChromosomeRef, Fixture
from .orientation_space import OrientationBasis, orientation_basis
from .order_dp import optimize_species_component_order


@dataclass(frozen=True)
class LocalSearchDiagnostics:
    restarts: int
    improving_steps: int
    candidate_evaluations: int
    order_dp_states_evaluated: int
    seed: int

    def to_dict(self) -> dict[str, object]:
        return {
            "restarts": self.restarts,
            "improving_steps": self.improving_steps,
            "candidate_evaluations": self.candidate_evaluations,
            "order_dp_states_evaluated": self.order_dp_states_evaluated,
            "seed": self.seed,
        }


@dataclass(frozen=True)
class LocalSearchResult:
    layout: ExactLayoutResult
    diagnostics: LocalSearchDiagnostics

    def to_dict(self) -> dict[str, object]:
        output = self.layout.to_dict()
        output["solver"] = "constraint-aware-local-search"
        output["diagnostics"] = self.diagnostics.to_dict()
        return output


@dataclass(frozen=True)
class AutoLayoutResult:
    layout: ExactLayoutResult
    solver: str
    details: dict[str, object]

    def to_dict(self) -> dict[str, object]:
        output = self.layout.to_dict()
        output["solver"] = self.solver
        output.update(self.details)
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


def _state_key(fixture: Fixture, state: LayoutState) -> tuple[object, ...]:
    return (
        tuple(
            ref.label
            for species in fixture.species_ids
            for ref in state.chromosome_order[species]
        ),
        tuple(
            state.chromosome_orientation[ref]
            for ref in sorted(state.chromosome_orientation)
        ),
    )


def _random_legal_state(
    fixture: Fixture,
    normalized: LayoutState,
    basis: OrientationBasis,
    component_of: dict[ChromosomeRef, int],
    rng: random.Random,
    randomize: bool,
) -> LayoutState:
    order: dict[str, tuple[ChromosomeRef, ...]] = {}

    for species in fixture.species_ids:
        current = normalized.chromosome_order[species]
        by_component: dict[int, list[ChromosomeRef]] = {}
        for ref in current:
            by_component.setdefault(component_of[ref], []).append(ref)

        rebuilt: list[ChromosomeRef] = []
        for component_id in sorted(by_component):
            refs = list(by_component[component_id])
            if randomize:
                rng.shuffle(refs)
            rebuilt.extend(refs)
        order[species] = tuple(rebuilt)

    orientation = dict(basis.base_assignment)
    if randomize:
        for group in basis.free_flip_groups:
            if rng.random() < 0.5:
                for ref in group:
                    orientation[ref] *= -1

    return LayoutState(
        chromosome_order=order,
        chromosome_orientation=orientation,
    )


def _candidate_moves(
    fixture: Fixture,
    state: LayoutState,
    component_of: dict[ChromosomeRef, int],
    basis: OrientationBasis,
):
    for species in fixture.species_ids:
        refs = state.chromosome_order[species]
        for index in range(len(refs) - 1):
            if component_of[refs[index]] == component_of[refs[index + 1]]:
                yield ("swap", species, index)

    for group_index, group in enumerate(basis.free_flip_groups):
        yield ("flip", group_index, group)


def _apply_move(state: LayoutState, move) -> LayoutState:
    kind = move[0]

    if kind == "swap":
        _, species, index = move
        order = dict(state.chromosome_order)
        refs = list(order[species])
        refs[index], refs[index + 1] = refs[index + 1], refs[index]
        order[species] = tuple(refs)
        return LayoutState(
            chromosome_order=order,
            chromosome_orientation=dict(state.chromosome_orientation),
        )

    _, _, group = move
    orientation = dict(state.chromosome_orientation)
    for ref in group:
        orientation[ref] *= -1
    return LayoutState(
        chromosome_order=dict(state.chromosome_order),
        chromosome_orientation=orientation,
    )


def optimize_local_search(
    fixture: Fixture,
    *,
    restarts: int = 8,
    max_improving_steps: int = 10000,
    seed: int = 1,
    progress_callback: Callable[[LayoutState, int], None] | None = None,
    order_dp_max_chromosomes: int | None = None,
    starting_state: LayoutState | None = None,
    prepared_components=None,
    prepared_basis: OrientationBasis | None = None,
) -> LocalSearchResult:
    """Constraint-aware best-improvement search over legal whole chromosomes.

    Each iteration can replace one species/component order with its exact
    O(n 2^n) subset-DP optimum conditional on neighboring layers. Adjacent
    whole-chromosome swaps remain available as a small local move. Orientation
    moves flip one complete GF(2) free group, preserving all hard equations.
    An optional chromosome cap skips expensive conditional DP moves only;
    skipped orders remain unresolved, not pruned from the legal search space.
    """

    if restarts < 1:
        raise ValueError("restarts must be at least 1")

    initial = initial_layout_state(fixture)
    normalized = canonicalize_component_order(fixture, initial)
    initial_score = score_crossings(fixture, initial)
    normalized_score = score_crossings(fixture, normalized)
    components, component_of = prepared_components or _component_map(fixture)
    basis = prepared_basis or orientation_basis(fixture, fixture.chromosome_refs)
    rng = random.Random(seed)

    best_state = None
    best_score = None
    evaluations = 0
    improving_steps = 0
    order_dp_states = 0
    cache = CrossingCostCache(fixture)
    reported_best = float("inf")

    def report(state, score):
        nonlocal reported_best
        if score < reported_best:
            reported_best = score
            if progress_callback is not None:
                progress_callback(state, score)

    for restart in range(restarts):
        state = _random_legal_state(
            fixture,
            normalized,
            basis,
            component_of,
            rng,
            randomize=(restart > 0),
        )
        if restart == 0 and starting_state is not None:
            # Validate external/warm starts through the canonical scorer.
            score_crossings(fixture, starting_state)
            for group in basis.free_flip_groups:
                signs = {starting_state.chromosome_orientation[ref] * basis.base_assignment[ref]
                         for ref in group}
                if len(signs) != 1:
                    raise ValueError("Warm start violates hard orientation equations")
            state = canonicalize_component_order(fixture, starting_state)
        score = cache.prepare(state)
        report(state, score)
        evaluations += 1
        steps = 0

        while steps < max_improving_steps:
            chosen = None

            for component_id, component_nodes in enumerate(components):
                for species in fixture.species_ids:
                    count = sum(component_of[ref] == component_id
                                for ref in state.chromosome_order[species])
                    if order_dp_max_chromosomes is not None and count > order_dp_max_chromosomes:
                        continue  # Work limit only; no order is declared excluded.
                    candidate, subproblem = optimize_species_component_order(
                        fixture, state, species, component_nodes
                    )
                    order_dp_states += subproblem.subset_states_evaluated
                    if candidate == state:
                        continue
                    candidate_score = cache.score_candidate(candidate)
                    evaluations += 1
                    move = ("reorder", species, component_id)
                    key = (
                        candidate_score,
                        move,
                        _state_key(fixture, candidate),
                    )
                    if candidate_score < score and (
                        chosen is None or key < chosen[0]
                    ):
                        chosen = (key, candidate, candidate_score)

            for move in _candidate_moves(fixture, state, component_of, basis):
                candidate = _apply_move(state, move)
                candidate_score = cache.score_candidate(candidate)
                evaluations += 1
                key = (candidate_score, move, _state_key(fixture, candidate))
                if candidate_score < score and (
                    chosen is None or key < chosen[0]
                ):
                    chosen = (key, candidate, candidate_score)

            if chosen is None:
                break

            state = chosen[1]
            score = chosen[2]
            if cache.prepare(state) != score:
                raise AssertionError("Cached candidate crossing score mismatch")
            report(state, score)
            steps += 1
            improving_steps += 1

        candidate_key = (score, _state_key(fixture, state))
        if best_score is None or candidate_key < best_score:
            best_score = candidate_key
            best_state = state

    assert best_state is not None
    optimized_score = score_crossings(fixture, best_state)
    if optimized_score.crossings != best_score[0]:
        raise AssertionError("Cached local-search score differs from canonical scorer")

    layout = ExactLayoutResult(
        initial_state=initial,
        normalized_state=normalized,
        optimized_state=best_state,
        initial_score=initial_score,
        normalized_score=normalized_score,
        optimized_score=optimized_score,
        states_evaluated=evaluations,
        optimality_status="best known",
    )
    return LocalSearchResult(
        layout=layout,
        diagnostics=LocalSearchDiagnostics(
            restarts=restarts,
            improving_steps=improving_steps,
            candidate_evaluations=evaluations,
            order_dp_states_evaluated=order_dp_states,
            seed=seed,
        ),
    )


def optimize_auto(
    fixture: Fixture,
    *,
    orientation_cap_per_component: int = 4096,
    permutation_cap_per_species: int = 40320,
    transition_cap_per_component: int = 5_000_000,
    branch_node_cap_per_component: int = 250_000,
    local_restarts: int = 8,
    local_max_improving_steps: int = 10000,
    seed: int = 1,
    component_workers: int = 1,
    time_limit_seconds: float | None = None,
    starting_state: LayoutState | None = None,
    bound_cluster_size: int = 0,
    progress_callback: Callable[[LayoutState, int], None] | None = None,
) -> AutoLayoutResult:
    """One component pipeline, with bounded continuation on reduced factors."""
    from .search_pipeline import optimize_pipeline
    return optimize_pipeline(
        fixture, orientation_cap_per_component=orientation_cap_per_component,
        permutation_cap_per_species=permutation_cap_per_species,
        transition_cap_per_component=transition_cap_per_component,
        branch_node_cap_per_component=branch_node_cap_per_component,
        local_restarts=local_restarts, local_max_improving_steps=local_max_improving_steps,
        seed=seed, component_workers=component_workers, progress_callback=progress_callback,
        time_limit_seconds=time_limit_seconds, starting_state=starting_state,
        bound_cluster_size=bound_cluster_size,
    )
