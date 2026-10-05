from __future__ import annotations

from dataclasses import dataclass
import random
from typing import Callable

from .crossing_cache import CrossingCostCache

from .incidence import build_incidence_graph, chromosome_node_id
from .layout import (
    ExactLayoutResult,
    LayoutState,
    SearchSpaceTooLarge,
    canonicalize_component_order,
    initial_layout_state,
    score_crossings,
)
from .layer_dp import LayerDPResult, exact_optimize_layer_dp
from .residual_solver import (
    ResidualExactResult,
    exact_optimize_residual_factor_graph,
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
) -> LocalSearchResult:
    """Constraint-aware best-improvement search over legal whole chromosomes.

    Each iteration can replace one species/component order with its exact
    O(n 2^n) subset-DP optimum conditional on neighboring layers. Adjacent
    whole-chromosome swaps remain available as a small local move. Orientation
    moves flip one complete GF(2) free group, preserving all hard equations.
    """

    if restarts < 1:
        raise ValueError("restarts must be at least 1")

    initial = initial_layout_state(fixture)
    normalized = canonicalize_component_order(fixture, initial)
    initial_score = score_crossings(fixture, initial)
    normalized_score = score_crossings(fixture, normalized)
    components, component_of = _component_map(fixture)
    basis = orientation_basis(fixture, fixture.chromosome_refs)
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
        score = cache.prepare(state)
        report(state, score)
        evaluations += 1
        steps = 0

        while steps < max_improving_steps:
            chosen = None

            for component_id, component_nodes in enumerate(components):
                for species in fixture.species_ids:
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
    progress_callback: Callable[[LayoutState, int], None] | None = None,
) -> AutoLayoutResult:
    """Use recursive residual elimination, then legacy exact DP/B&B fallback.

    The residual solver is now the first exact method because it repeatedly
    factors the *remaining objective* rather than enumerating global
    orientation assignments. If an intermediate factor would exceed the same
    configured transition/work cap, the historical exact layer DP is still
    attempted before bounded branch-and-bound.
    """

    residual_failure: str | None = None
    try:
        residual: ResidualExactResult = (
            exact_optimize_residual_factor_graph(
                fixture,
                permutation_cap_per_variable=permutation_cap_per_species,
                table_entry_cap_per_component=transition_cap_per_component,
                work_cap_per_component=transition_cap_per_component,
                component_workers=component_workers,
            )
        )
        return AutoLayoutResult(
            layout=residual.layout,
            solver="exact-residual-factor-elimination",
            details={
                "component_diagnostics": [
                    item.to_dict() for item in residual.diagnostics
                ],
                "component_workers": residual.component_workers,
            },
        )
    except SearchSpaceTooLarge as exc:
        residual_failure = str(exc)

    try:
        exact: LayerDPResult = exact_optimize_layer_dp(
            fixture,
            orientation_cap_per_component=orientation_cap_per_component,
            permutation_cap_per_species=permutation_cap_per_species,
            transition_cap_per_component=transition_cap_per_component,
        )
        return AutoLayoutResult(
            layout=exact.layout,
            solver="exact-layer-dynamic-programming",
            details={
                "residual_fallback_reason": residual_failure,
                "component_diagnostics": [
                    item.to_dict() for item in exact.diagnostics
                ],
            },
        )
    except SearchSpaceTooLarge as exc:
        # Local import avoids a module-level cycle: branch_bound uses the
        # local-search routine above as its incumbent generator.
        from .branch_bound import optimize_branch_and_bound

        bounded = optimize_branch_and_bound(
            fixture,
            node_cap_per_component=branch_node_cap_per_component,
            orientation_cap_per_component=orientation_cap_per_component,
            local_restarts=local_restarts,
            local_max_improving_steps=local_max_improving_steps,
            progress_callback=progress_callback,
            seed=seed,
            component_workers=component_workers,
        )
        return AutoLayoutResult(
            layout=bounded.layout,
            solver="monotone-component-branch-and-bound",
            details={
                "residual_fallback_reason": residual_failure,
                "fallback_reason": str(exc),
                "lower_bound": bounded.lower_bound,
                "upper_bound": bounded.upper_bound,
                "optimality_gap": bounded.gap,
                "component_diagnostics": [
                    item.to_dict() for item in bounded.components
                ],
                "component_workers": bounded.component_workers,
            },
        )
