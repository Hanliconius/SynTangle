from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
import multiprocessing as mp
from itertools import permutations, product
from math import factorial

from .incidence import build_incidence_graph, chromosome_node_id
from .layout import (
    ExactLayoutResult,
    LayoutState,
    SearchSpaceTooLarge,
    _occurrences_by_species_homology,
    canonicalize_component_order,
    initial_layout_state,
    score_crossings,
)
from .model import ChromosomeRef, Fixture
from .orientation_space import OrientationBasis, orientation_basis
from .pair_cost import pair_component_crossings
from .residual import (
    ResidualFactor,
    ResidualVariable,
    build_residual_factorization,
)


@dataclass(frozen=True)
class ResidualSolveDiagnostics:
    component_id: int
    chromosome_count: int
    variable_count: int
    factor_count: int
    objective_component_count: int
    isolated_variable_count: int
    leaf_eliminations: int
    articulation_conditionings: int
    min_fill_eliminations: int
    dynamic_factor_splits: int
    factor_scope_variables_removed: int
    max_intermediate_scope: int
    max_table_entries: int
    table_entries_evaluated: int
    component_workers: int = 1

    def to_dict(self) -> dict[str, object]:
        return {
            "component_id": self.component_id,
            "chromosome_count": self.chromosome_count,
            "variable_count": self.variable_count,
            "factor_count": self.factor_count,
            "objective_component_count": self.objective_component_count,
            "isolated_variable_count": self.isolated_variable_count,
            "leaf_eliminations": self.leaf_eliminations,
            "articulation_conditionings": self.articulation_conditionings,
            "min_fill_eliminations": self.min_fill_eliminations,
            "dynamic_factor_splits": self.dynamic_factor_splits,
            "factor_scope_variables_removed": (
                self.factor_scope_variables_removed
            ),
            "max_intermediate_scope": self.max_intermediate_scope,
            "max_table_entries": self.max_table_entries,
            "table_entries_evaluated": self.table_entries_evaluated,
            "component_workers": self.component_workers,
        }


@dataclass(frozen=True)
class ResidualExactResult:
    layout: ExactLayoutResult
    diagnostics: tuple[ResidualSolveDiagnostics, ...]
    component_workers: int = 1

    def to_dict(self) -> dict[str, object]:
        output = self.layout.to_dict()
        output["solver"] = "exact-residual-factor-elimination"
        output["component_diagnostics"] = [
            item.to_dict() for item in self.diagnostics
        ]
        output["component_workers"] = self.component_workers
        return output


@dataclass(frozen=True)
class _Domain:
    variable_id: str
    values: tuple[object, ...]

    @property
    def size(self) -> int:
        return len(self.values)


@dataclass(frozen=True)
class _TableFactor:
    factor_id: str
    scope: tuple[str, ...]
    values: dict[tuple[int, ...], int]


@dataclass(frozen=True)
class _EliminationRecord:
    variable_id: str
    remaining_scope: tuple[str, ...]
    best_value: dict[tuple[int, ...], int]


@dataclass
class _WorkBudget:
    cap: int
    table_cap: int
    used: int = 0
    max_table_entries: int = 0
    max_intermediate_scope: int = 0
    continuation_node_cap: int | None = None
    continuation_nodes: int = 0
    continuation_pruned: int = 0
    preferred: dict[str, int] = field(default_factory=dict)
    handoffs: list[dict[str, object]] = field(default_factory=list)
    reconstruction_stack: list[_EliminationRecord] = field(default_factory=list)
    conditions: dict[str, int] = field(default_factory=dict)
    scope_reductions: list[dict[str, object]] = field(default_factory=list)

    def consume(self, count: int, *, scope_size: int) -> None:
        if count < 0:
            raise ValueError("Work count cannot be negative")
        if count > self.table_cap:
            raise SearchSpaceTooLarge(
                "Residual factor elimination would create/evaluate a table "
                f"with {count} assignments; table cap is {self.table_cap}"
            )
        if self.used + count > self.cap:
            raise SearchSpaceTooLarge(
                "Residual factor elimination would exceed the configured "
                f"work cap ({self.used + count} > {self.cap})"
            )
        self.used += count
        self.max_table_entries = max(self.max_table_entries, count)
        self.max_intermediate_scope = max(
            self.max_intermediate_scope,
            scope_size,
        )


@dataclass
class _MutableDiagnostics:
    leaf_eliminations: int = 0
    articulation_conditionings: int = 0
    min_fill_eliminations: int = 0
    dynamic_factor_splits: int = 0
    factor_scope_variables_removed: int = 0


@dataclass(frozen=True)
class _FactorSolveResult:
    cost: int
    assignment: dict[str, int]
    lower_bound: int | None = None

    @property
    def lower(self) -> int:
        return self.cost if self.lower_bound is None else self.lower_bound


@dataclass(frozen=True)
class _ComponentResidualResult:
    component_id: int
    orders: dict[str, tuple[ChromosomeRef, ...]]
    orientation: dict[ChromosomeRef, int]
    optimum: int
    diagnostics: ResidualSolveDiagnostics
    lower_bound: int = 0
    handoffs: tuple[dict[str, object], ...] = ()
    continuation_nodes: int = 0
    continuation_pruned: int = 0


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


def _assignment_key(assignment: dict[str, int]) -> tuple[tuple[str, int], ...]:
    return tuple(sorted(assignment.items()))


def _table_size(
    scope: tuple[str, ...],
    domains: dict[str, _Domain],
) -> int:
    size = 1
    for variable_id in scope:
        size *= domains[variable_id].size
    return size


def _factor_value(
    factor: _TableFactor,
    assignment: dict[str, int],
) -> int:
    return factor.values[
        tuple(assignment[variable_id] for variable_id in factor.scope)
    ]


def _extract_constants(
    factors: list[_TableFactor],
) -> tuple[int, list[_TableFactor]]:
    constant = 0
    remaining: list[_TableFactor] = []
    for factor in factors:
        if factor.scope:
            remaining.append(factor)
        else:
            constant += factor.values[()]
    return constant, remaining


def _factor_components(
    factors: list[_TableFactor],
) -> list[list[_TableFactor]]:
    active = [factor for factor in factors if factor.scope]
    if not active:
        return []

    by_variable: dict[str, list[int]] = {}
    for index, factor in enumerate(active):
        for variable_id in factor.scope:
            by_variable.setdefault(variable_id, []).append(index)

    unseen = set(range(len(active)))
    components: list[list[_TableFactor]] = []

    while unseen:
        seed = min(unseen)
        stack = [seed]
        seen_factors: set[int] = set()
        seen_variables: set[str] = set()

        while stack:
            index = stack.pop()
            if index in seen_factors:
                continue
            seen_factors.add(index)
            unseen.discard(index)
            factor = active[index]
            for variable_id in factor.scope:
                if variable_id in seen_variables:
                    continue
                seen_variables.add(variable_id)
                for neighbor in by_variable[variable_id]:
                    if neighbor not in seen_factors:
                        stack.append(neighbor)

        components.append([active[index] for index in sorted(seen_factors)])

    components.sort(
        key=lambda group: (
            -len({v for factor in group for v in factor.scope}),
            tuple(factor.factor_id for factor in group),
        )
    )
    return components


def _primal_adjacency(
    factors: list[_TableFactor],
) -> dict[str, set[str]]:
    variables = {
        variable_id
        for factor in factors
        for variable_id in factor.scope
    }
    adjacency = {variable_id: set() for variable_id in variables}

    for factor in factors:
        for index, a in enumerate(factor.scope):
            for b in factor.scope[index + 1 :]:
                adjacency[a].add(b)
                adjacency[b].add(a)
    return adjacency


def _articulation_points(
    adjacency: dict[str, set[str]],
) -> tuple[str, ...]:
    discovery: dict[str, int] = {}
    low: dict[str, int] = {}
    parent: dict[str, str | None] = {}
    points: set[str] = set()
    counter = 0

    def visit(node: str) -> None:
        nonlocal counter
        counter += 1
        discovery[node] = counter
        low[node] = counter
        children = 0

        for neighbor in sorted(adjacency[node]):
            if neighbor not in discovery:
                parent[neighbor] = node
                children += 1
                visit(neighbor)
                low[node] = min(low[node], low[neighbor])
                if parent.get(node) is None and children > 1:
                    points.add(node)
                if (
                    parent.get(node) is not None
                    and low[neighbor] >= discovery[node]
                ):
                    points.add(node)
            elif neighbor != parent.get(node):
                low[node] = min(low[node], discovery[neighbor])

    for node in sorted(adjacency):
        if node in discovery:
            continue
        parent[node] = None
        visit(node)

    return tuple(sorted(points))


def _condition_factor(
    factor: _TableFactor,
    variable_id: str,
    value: int,
) -> _TableFactor:
    if variable_id not in factor.scope:
        return factor

    index = factor.scope.index(variable_id)
    scope = tuple(
        item for item in factor.scope if item != variable_id
    )
    values: dict[tuple[int, ...], int] = {}

    for key, cost in factor.values.items():
        if key[index] != value:
            continue
        reduced = key[:index] + key[index + 1 :]
        values[reduced] = cost

    return _TableFactor(
        factor_id=f"{factor.factor_id}|{variable_id}={value}",
        scope=scope,
        values=values,
    )


def _eliminate_variable(
    variable_id: str,
    factors: list[_TableFactor],
    domains: dict[str, _Domain],
    budget: _WorkBudget,
) -> tuple[list[_TableFactor], _EliminationRecord]:
    involved = [
        factor for factor in factors if variable_id in factor.scope
    ]
    if not involved:
        return (
            factors,
            _EliminationRecord(
                variable_id=variable_id,
                remaining_scope=(),
                best_value={(): 0},
            ),
        )

    untouched = [
        factor for factor in factors if variable_id not in factor.scope
    ]
    union_scope = tuple(
        sorted(
            {
                item
                for factor in involved
                for item in factor.scope
            }
        )
    )
    remaining_scope = tuple(
        item for item in union_scope if item != variable_id
    )

    entry_count = _table_size(union_scope, domains)
    budget.consume(entry_count, scope_size=len(union_scope))

    ranges = [
        range(domains[item].size)
        for item in union_scope
    ]
    best_cost: dict[tuple[int, ...], int] = {}
    best_value: dict[tuple[int, ...], int] = {}

    for values in product(*ranges):
        assignment = dict(zip(union_scope, values))
        cost = sum(
            _factor_value(factor, assignment)
            for factor in involved
        )
        remaining_key = tuple(
            assignment[item] for item in remaining_scope
        )
        candidate_value = assignment[variable_id]

        previous = best_cost.get(remaining_key)
        if (
            previous is None
            or cost < previous
            or (
                cost == previous
                and candidate_value < best_value[remaining_key]
            )
        ):
            best_cost[remaining_key] = cost
            best_value[remaining_key] = candidate_value

    reduced = _TableFactor(
        factor_id=f"elim::{variable_id}",
        scope=remaining_scope,
        values=best_cost,
    )
    return (
        untouched + [reduced],
        _EliminationRecord(
            variable_id=variable_id,
            remaining_scope=remaining_scope,
            best_value=best_value,
        ),
    )


def _min_fill_variable(
    factors: list[_TableFactor],
    domains: dict[str, _Domain],
) -> str:
    adjacency = _primal_adjacency(factors)

    def key(variable_id: str) -> tuple[int, int, int, str]:
        neighbors = sorted(adjacency[variable_id])
        fill = 0
        for index, a in enumerate(neighbors):
            for b in neighbors[index + 1 :]:
                if b not in adjacency[a]:
                    fill += 1

        induced_entries = domains[variable_id].size
        for neighbor in neighbors:
            induced_entries *= domains[neighbor].size

        return (
            fill,
            induced_entries,
            len(neighbors),
            variable_id,
        )

    return min(adjacency, key=key)


def _bounded_factor_search(factors, domains, budget, diagnostics):
    """Continue on precisely the current factors; never rebuild original scopes.

    Eliminations above this call remain on the reconstruction stack. Conditioning
    is already encoded in these tables, and disconnected solved pieces remain
    with their caller. A node cap leaves a valid bound, never an exclusion.
    """
    import heapq
    variables = tuple(sorted({v for f in factors for v in f.scope}))
    eliminated = {r.variable_id for r in budget.reconstruction_stack}
    if eliminated.intersection(variables) or set(budget.conditions).intersection(variables):
        raise AssertionError("Eliminated/conditioned variable reintroduced at handoff")
    budget.handoffs.append({
        "active_variables": list(variables),
        "retained_eliminations": [dict(variable=r.variable_id,
            remaining_scope=list(r.remaining_scope), reason="exact conditional minimization",
            reconstruction_entries=len(r.best_value)) for r in budget.reconstruction_stack],
        "retained_conditions": dict(budget.conditions),
        "retained_scope_reductions": list(budget.scope_reductions),
        "factor_scopes": [list(f.scope) for f in factors],
        "leaf_eliminations_retained": diagnostics.leaf_eliminations,
        "min_fill_eliminations_retained": diagnostics.min_fill_eliminations,
        "scope_variables_removed_retained": diagnostics.factor_scope_variables_removed,
        "reason": "elimination work/table cap; continue reduced factors",
        "scope": "current conditioned residual piece",
    })
    preferred = {v: budget.preferred.get(v, 0) for v in variables}
    best = dict(preferred)
    upper = sum(_factor_value(f, best) for f in factors)

    def lower(partial):
        total = 0
        for factor in factors:
            positions = [(i, partial[v]) for i, v in enumerate(factor.scope) if v in partial]
            total += min(cost for key, cost in factor.values.items()
                         if all(key[i] == value for i, value in positions))
        return total

    root_lower = lower({})
    serial = 0
    heap = [(root_lower, serial, {})]
    while heap and budget.continuation_nodes < budget.continuation_node_cap:
        bound, _, partial = heapq.heappop(heap)
        if bound >= upper:
            budget.continuation_pruned += 1
            continue
        budget.continuation_nodes += 1
        if len(partial) == len(variables):
            upper, best = bound, partial
            continue
        variable = min((v for v in variables if v not in partial),
                       key=lambda v: (domains[v].size, v))
        values = [preferred[variable]] + [i for i in range(domains[variable].size)
                                         if i != preferred[variable]]
        for value in values:
            child = {**partial, variable: value}
            child_lower = lower(child)
            if child_lower >= upper:
                budget.continuation_pruned += 1
                continue
            feasible = {**preferred, **child}
            cost = sum(_factor_value(f, feasible) for f in factors)
            if cost < upper:
                upper, best = cost, feasible
            if child_lower < upper:
                serial += 1
                heapq.heappush(heap, (child_lower, serial, child))
            else:
                budget.continuation_pruned += 1
    final_lower = min(upper, heap[0][0]) if heap else upper
    budget.handoffs[-1].update(lower_bound=final_lower, upper_bound=upper,
                              unresolved=final_lower < upper)
    return _FactorSolveResult(upper, best, final_lower)


def _solve_factor_system(
    factors: list[_TableFactor],
    domains: dict[str, _Domain],
    budget: _WorkBudget,
    diagnostics: _MutableDiagnostics,
    *,
    separator_domain_cap: int,
) -> _FactorSolveResult:
    constant, active_factors = _extract_constants(factors)

    if not active_factors:
        return _FactorSolveResult(cost=constant, assignment={})

    pieces = _factor_components(active_factors)
    if len(pieces) > 1:
        diagnostics.dynamic_factor_splits += len(pieces) - 1
        total = constant
        total_lower = constant
        assignment: dict[str, int] = {}
        for piece in pieces:
            result = _solve_factor_system(
                piece,
                domains,
                budget,
                diagnostics,
                separator_domain_cap=separator_domain_cap,
            )
            total += result.cost
            total_lower += result.lower
            overlap = set(assignment) & set(result.assignment)
            if overlap:
                raise AssertionError(
                    "Disconnected residual pieces unexpectedly share variables"
                )
            assignment.update(result.assignment)
        return _FactorSolveResult(total, assignment, total_lower)

    active_factors = pieces[0]
    adjacency = _primal_adjacency(active_factors)

    # "Leaf" is defined on the residual primal graph, not by the raw number of
    # factor records containing a variable. After eliminating one endpoint of
    # a chain, its neighbor commonly has both a unary message factor and one
    # pair factor: two factor occurrences but only one remaining neighbor. Such
    # a variable is still an exact leaf and should be peeled rather than
    # triggering separator branching.
    leaf_variables = sorted(
        variable_id
        for variable_id, neighbors in adjacency.items()
        if len(neighbors) <= 1
    )
    if leaf_variables:
        variable_id = min(
            leaf_variables,
            key=lambda item: (
                _table_size(
                    tuple(
                        sorted(
                            {
                                v
                                for factor in active_factors
                                if item in factor.scope
                                for v in factor.scope
                            }
                        )
                    ),
                    domains,
                ),
                domains[item].size,
                item,
            ),
        )
        try:
            reduced, record = _eliminate_variable(
                variable_id, active_factors, domains, budget,
            )
        except SearchSpaceTooLarge:
            if budget.continuation_node_cap is None:
                raise
            result = _bounded_factor_search(active_factors, domains, budget, diagnostics)
            return _FactorSolveResult(constant + result.cost, result.assignment,
                                      constant + result.lower)
        diagnostics.leaf_eliminations += 1
        budget.reconstruction_stack.append(record)
        try:
            result = _solve_factor_system(
                reduced, domains, budget, diagnostics,
                separator_domain_cap=separator_domain_cap,
            )
        finally:
            budget.reconstruction_stack.pop()
        remaining_key = tuple(
            result.assignment[item]
            for item in record.remaining_scope
        )
        assignment = dict(result.assignment)
        assignment[variable_id] = record.best_value[remaining_key]
        return _FactorSolveResult(
            cost=constant + result.cost,
            assignment=assignment,
            lower_bound=constant + result.lower,
        )

    articulation = [
        variable_id
        for variable_id in _articulation_points(adjacency)
        if domains[variable_id].size <= separator_domain_cap
    ]
    if articulation:
        variable_id = min(
            articulation,
            key=lambda item: (
                domains[item].size,
                len(adjacency[item]),
                item,
            ),
        )
        diagnostics.articulation_conditionings += 1
        best: tuple[
            tuple[int, tuple[tuple[str, int], ...]],
            _FactorSolveResult,
        ] | None = None

        branch_lowers = []
        for value in range(domains[variable_id].size):
            conditioned = [
                _condition_factor(factor, variable_id, value)
                for factor in active_factors
            ]
            budget.conditions[variable_id] = value
            try:
                result = _solve_factor_system(
                    conditioned, domains, budget, diagnostics,
                    separator_domain_cap=separator_domain_cap,
                )
            finally:
                del budget.conditions[variable_id]
            branch_lowers.append(constant + result.lower)
            assignment = dict(result.assignment)
            assignment[variable_id] = value
            candidate = _FactorSolveResult(
                cost=constant + result.cost,
                assignment=assignment,
            )
            key = (candidate.cost, _assignment_key(candidate.assignment))
            if best is None or key < best[0]:
                best = (key, candidate)

        assert best is not None
        return _FactorSolveResult(best[1].cost, best[1].assignment, min(branch_lowers))

    variable_id = _min_fill_variable(active_factors, domains)
    try:
        reduced, record = _eliminate_variable(variable_id, active_factors, domains, budget)
    except SearchSpaceTooLarge:
        if budget.continuation_node_cap is None:
            raise
        result = _bounded_factor_search(active_factors, domains, budget, diagnostics)
        return _FactorSolveResult(constant + result.cost, result.assignment,
                                  constant + result.lower)
    diagnostics.min_fill_eliminations += 1
    budget.reconstruction_stack.append(record)
    try:
        result = _solve_factor_system(
            reduced, domains, budget, diagnostics,
            separator_domain_cap=separator_domain_cap,
        )
    finally:
        budget.reconstruction_stack.pop()
    remaining_key = tuple(
        result.assignment[item]
        for item in record.remaining_scope
    )
    assignment = dict(result.assignment)
    assignment[variable_id] = record.best_value[remaining_key]
    return _FactorSolveResult(
        cost=constant + result.cost,
        assignment=assignment,
        lower_bound=constant + result.lower,
    )


def _orientation_from_assignment(
    basis: OrientationBasis,
    orientation_variables: tuple[ResidualVariable, ...],
    assignment: dict[str, int],
) -> dict[ChromosomeRef, int]:
    orientation = dict(basis.base_assignment)
    for variable in orientation_variables:
        if assignment.get(variable.variable_id, 0) == 0:
            continue
        for ref in variable.chromosome_refs:
            orientation[ref] *= -1
    return orientation


def _component_orders(
    fixture: Fixture,
    refs: tuple[ChromosomeRef, ...],
    order_variables: tuple[ResidualVariable, ...],
    domains: dict[str, _Domain],
    assignment: dict[str, int],
) -> dict[str, tuple[ChromosomeRef, ...]]:
    variable_by_species = {
        variable.species_id: variable
        for variable in order_variables
    }
    orders: dict[str, tuple[ChromosomeRef, ...]] = {}

    for species in fixture.species_ids:
        species_refs = tuple(
            ref for ref in refs if ref.species_id == species
        )
        variable = variable_by_species.get(species)
        if variable is None:
            orders[species] = species_refs
            continue

        value_index = assignment.get(variable.variable_id, 0)
        value = domains[variable.variable_id].values[value_index]
        if not isinstance(value, tuple):
            raise AssertionError("Order domain contains a non-tuple value")
        orders[species] = value

    return orders


def _build_domains(
    variables: tuple[ResidualVariable, ...],
    *,
    permutation_cap_per_variable: int,
) -> dict[str, _Domain]:
    domains: dict[str, _Domain] = {}

    for variable in variables:
        if variable.kind == "orientation":
            values: tuple[object, ...] = (0, 1)
        elif variable.kind == "order":
            count = factorial(len(variable.chromosome_refs))
            if count > permutation_cap_per_variable:
                raise SearchSpaceTooLarge(
                    f"Residual order variable {variable.variable_id} has "
                    f"{count} permutations; cap is "
                    f"{permutation_cap_per_variable}"
                )
            values = tuple(permutations(variable.chromosome_refs))
        else:
            raise ValueError(
                f"Unknown residual variable kind: {variable.kind}"
            )

        domains[variable.variable_id] = _Domain(
            variable_id=variable.variable_id,
            values=values,
        )

    return domains


def _evaluate_residual_factor(
    fixture: Fixture,
    factor: ResidualFactor,
    component_nodes: frozenset[str],
    refs: tuple[ChromosomeRef, ...],
    basis: OrientationBasis,
    variables_by_id: dict[str, ResidualVariable],
    domains: dict[str, _Domain],
    assignment: dict[str, int],
    occurrence_index,
) -> int:
    order_variables = tuple(
        variable
        for variable in variables_by_id.values()
        if variable.kind == "order"
    )
    orientation_variables = tuple(
        variable
        for variable in variables_by_id.values()
        if variable.kind == "orientation"
    )

    orders = _component_orders(
        fixture,
        refs,
        order_variables,
        domains,
        assignment,
    )
    orientation = _orientation_from_assignment(
        basis,
        orientation_variables,
        assignment,
    )
    return pair_component_crossings(
        fixture,
        factor.species_left,
        factor.species_right,
        orders[factor.species_left],
        orders[factor.species_right],
        orientation,
        component_nodes,
        occurrence_index=occurrence_index,
    )


def _drop_invariant_variables(
    factor: _TableFactor,
    domains: dict[str, _Domain],
) -> tuple[_TableFactor, int]:
    """Remove conservative scope variables that provably do not affect a table.

    Residual graph construction deliberately over-approximates dependencies.
    Once a factor table has been evaluated, exact value equality lets us
    contract any variable whose value leaves that factor unchanged. Repeating
    this to a fixed point can expose new disconnected pieces and leaves.
    """

    scope = list(factor.scope)
    values = dict(factor.values)
    removed = 0

    changed = True
    while changed and scope:
        changed = False
        for index, variable_id in enumerate(tuple(scope)):
            other_scope = tuple(
                item for item in scope if item != variable_id
            )
            grouped: dict[tuple[int, ...], list[int]] = {}

            for key, cost in values.items():
                reduced_key = key[:index] + key[index + 1 :]
                grouped.setdefault(reduced_key, []).append(cost)

            if any(
                len(costs) != domains[variable_id].size
                or len(set(costs)) != 1
                for costs in grouped.values()
            ):
                continue

            scope.pop(index)
            values = {
                reduced_key: costs[0]
                for reduced_key, costs in grouped.items()
            }
            removed += 1
            changed = True
            break

    return (
        _TableFactor(
            factor_id=factor.factor_id,
            scope=tuple(scope),
            values=values,
        ),
        removed,
    )


def _build_factor_tables(
    fixture: Fixture,
    component_nodes: frozenset[str],
    refs: tuple[ChromosomeRef, ...],
    basis: OrientationBasis,
    variables: tuple[ResidualVariable, ...],
    factors: tuple[ResidualFactor, ...],
    domains: dict[str, _Domain],
    budget: _WorkBudget,
    *,
    table_entry_cap: int,
) -> list[_TableFactor]:
    variables_by_id = {
        variable.variable_id: variable
        for variable in variables
    }
    occurrence_index = _occurrences_by_species_homology(fixture)
    output: list[_TableFactor] = []

    sizes = [_table_size(tuple(sorted(f.variable_ids)), domains) for f in factors]
    if any(size > table_entry_cap for size in sizes) or budget.used + sum(sizes) > budget.cap:
        raise SearchSpaceTooLarge("Factor materialization exceeds cap; no factor reductions started")

    for factor in factors:
        scope = tuple(sorted(factor.variable_ids))
        entry_count = _table_size(scope, domains)
        if entry_count > table_entry_cap:
            raise SearchSpaceTooLarge(
                f"Residual factor {factor.factor_id} needs "
                f"{entry_count} entries; cap is {table_entry_cap}"
            )
        budget.consume(entry_count, scope_size=len(scope))

        values: dict[tuple[int, ...], int] = {}
        ranges = [range(domains[item].size) for item in scope]
        iterator = product(*ranges) if ranges else [()]

        for key in iterator:
            assignment = dict(zip(scope, key))
            values[tuple(key)] = _evaluate_residual_factor(
                fixture,
                factor,
                component_nodes,
                refs,
                basis,
                variables_by_id,
                domains,
                assignment,
                occurrence_index,
            )

        table = _TableFactor(
            factor_id=factor.factor_id,
            scope=scope,
            values=values,
        )
        table, removed = _drop_invariant_variables(
            table, domains,
        )
        if removed:
            budget.scope_reductions.append(dict(
                factor=factor.factor_id, removed_variables=sorted(set(scope)-set(table.scope)),
                reason="exact table invariance", scope="this factor only"))
        output.append(table)

    return output


def _solve_incidence_component(
    fixture: Fixture,
    component_id: int,
    component_nodes: frozenset[str],
    component_of: dict[ChromosomeRef, int],
    variables: tuple[ResidualVariable, ...],
    factors: tuple[ResidualFactor, ...],
    permutation_cap_per_variable: int,
    table_entry_cap_per_component: int,
    work_cap_per_component: int,
    separator_domain_cap: int,
    *,
    prepared_basis: OrientationBasis | None = None,
    incumbent_state: LayoutState | None = None,
    continuation_node_cap: int | None = None,
) -> _ComponentResidualResult:
    refs = tuple(
        sorted(
            ref
            for ref in fixture.chromosome_refs
            if component_of[ref] == component_id
        )
    )
    basis = prepared_basis or orientation_basis(fixture, refs)
    domains = _build_domains(
        variables,
        permutation_cap_per_variable=permutation_cap_per_variable,
    )
    budget = _WorkBudget(
        cap=work_cap_per_component,
        table_cap=table_entry_cap_per_component,
        continuation_node_cap=continuation_node_cap,
    )
    if incumbent_state is not None:
        for variable in variables:
            if variable.kind == "orientation":
                ref = variable.chromosome_refs[0]
                budget.preferred[variable.variable_id] = int(
                    incumbent_state.chromosome_orientation[ref] != basis.base_assignment[ref])
            else:
                selected = tuple(ref for ref in incumbent_state.chromosome_order[variable.species_id]
                                 if ref in variable.chromosome_refs)
                budget.preferred[variable.variable_id] = domains[variable.variable_id].values.index(selected)
    mutable = _MutableDiagnostics()

    raw_tables = _build_factor_tables(
        fixture,
        component_nodes,
        refs,
        basis,
        variables,
        factors,
        domains,
        budget,
        table_entry_cap=table_entry_cap_per_component,
    )
    original_scope_size = sum(
        len(factor.variable_ids) for factor in factors
    )
    reduced_scope_size = sum(
        len(table.scope) for table in raw_tables
    )
    mutable.factor_scope_variables_removed = (
        original_scope_size - reduced_scope_size
    )
    tables = raw_tables

    active_variable_ids = {
        variable_id
        for factor in tables
        for variable_id in factor.scope
    }
    isolated = tuple(
        sorted(set(domains) - active_variable_ids)
    )

    initial_pieces = _factor_components(tables)
    result = _solve_factor_system(
        tables,
        domains,
        budget,
        mutable,
        separator_domain_cap=separator_domain_cap,
    )

    assignment = dict(result.assignment)
    for variable_id in isolated:
        assignment.setdefault(variable_id, 0)

    order_variables = tuple(
        variable for variable in variables if variable.kind == "order"
    )
    orientation_variables = tuple(
        variable
        for variable in variables
        if variable.kind == "orientation"
    )
    orders = _component_orders(
        fixture,
        refs,
        order_variables,
        domains,
        assignment,
    )
    orientation = _orientation_from_assignment(
        basis,
        orientation_variables,
        assignment,
    )

    # Verify that the factor-table optimum exactly matches direct scoring for
    # the materialized legal component state.
    full_order = {
        species: tuple(
            ref
            for ref in fixture.chromosome_refs
            if ref.species_id == species
            and component_of[ref] == component_id
        )
        for species in fixture.species_ids
    }
    full_order.update(orders)
    state = LayoutState(
        chromosome_order=full_order,
        chromosome_orientation=orientation,
    )
    direct = score_crossings(
        fixture,
        state,
        restrict_component_nodes=component_nodes,
    ).crossings
    if direct != result.cost:
        raise AssertionError(
            "Residual factor optimum does not match direct component "
            f"crossing score: factor={result.cost}, direct={direct}"
        )

    diagnostic = ResidualSolveDiagnostics(
        component_id=component_id,
        chromosome_count=len(refs),
        variable_count=len(variables),
        factor_count=len(factors),
        objective_component_count=len(initial_pieces),
        isolated_variable_count=len(isolated),
        leaf_eliminations=mutable.leaf_eliminations,
        articulation_conditionings=(
            mutable.articulation_conditionings
        ),
        min_fill_eliminations=mutable.min_fill_eliminations,
        dynamic_factor_splits=mutable.dynamic_factor_splits,
        factor_scope_variables_removed=(
            mutable.factor_scope_variables_removed
        ),
        max_intermediate_scope=budget.max_intermediate_scope,
        max_table_entries=budget.max_table_entries,
        table_entries_evaluated=budget.used,
    )
    return _ComponentResidualResult(
        component_id=component_id,
        orders=orders,
        orientation=orientation,
        optimum=result.cost,
        diagnostics=diagnostic,
        lower_bound=result.lower,
        handoffs=tuple(budget.handoffs),
        continuation_nodes=budget.continuation_nodes,
        continuation_pruned=budget.continuation_pruned,
    )


def exact_optimize_residual_factor_graph(
    fixture: Fixture,
    *,
    permutation_cap_per_variable: int = 40320,
    table_entry_cap_per_component: int = 250_000,
    work_cap_per_component: int = 5_000_000,
    separator_domain_cap: int = 8,
    component_workers: int = 1,
) -> ResidualExactResult:
    """Exactly solve the crossing objective by recursive residual reduction.

    The solver works on the decision/factor graph rather than enumerating one
    global orientation state. At every recursive step it:

    1. removes constant/objective-neutral structure;
    2. splits disconnected factor pieces exactly;
    3. eliminates variables occurring in a single factor;
    4. conditions on small articulation variables when that exposes
       independent subproblems; and
    5. otherwise performs a min-fill variable elimination step.

    Every transformation is exact. If an intermediate factor table would
    exceed the configured caps, SearchSpaceTooLarge is raised and callers may
    fall back to another exact/bounded solver.
    """

    if permutation_cap_per_variable < 1:
        raise ValueError("permutation_cap_per_variable must be at least 1")
    if table_entry_cap_per_component < 1:
        raise ValueError(
            "table_entry_cap_per_component must be at least 1"
        )
    if work_cap_per_component < 1:
        raise ValueError("work_cap_per_component must be at least 1")
    if separator_domain_cap < 1:
        raise ValueError("separator_domain_cap must be at least 1")
    if component_workers < 1:
        raise ValueError("component_workers must be at least 1")

    initial = initial_layout_state(fixture)
    normalized = canonicalize_component_order(fixture, initial)
    initial_score = score_crossings(fixture, initial)
    normalized_score = score_crossings(fixture, normalized)

    components, component_of = _component_map(fixture)
    residual = build_residual_factorization(fixture)

    variables_by_component = {
        component_id: tuple(
            variable
            for variable in residual.variables
            if variable.component_id == component_id
        )
        for component_id in range(len(components))
    }
    factors_by_component = {
        component_id: tuple(
            factor
            for factor in residual.factors
            if factor.component_id == component_id
        )
        for component_id in range(len(components))
    }

    worker_count = min(component_workers, max(1, len(components)))
    arguments = [
        (
            fixture,
            component_id,
            component_nodes,
            component_of,
            variables_by_component[component_id],
            factors_by_component[component_id],
            permutation_cap_per_variable,
            table_entry_cap_per_component,
            work_cap_per_component,
            separator_domain_cap,
        )
        for component_id, component_nodes in enumerate(components)
    ]

    if worker_count > 1 and len(arguments) > 1:
        with ProcessPoolExecutor(
            max_workers=worker_count,
            mp_context=mp.get_context("spawn"),
        ) as executor:
            futures = [
                executor.submit(_solve_incidence_component, *args)
                for args in arguments
            ]
            component_results = [
                future.result() for future in futures
            ]
    else:
        component_results = [
            _solve_incidence_component(*args)
            for args in arguments
        ]

    component_results.sort(key=lambda item: item.component_id)

    final_order: dict[str, tuple[ChromosomeRef, ...]] = {}
    for species in fixture.species_ids:
        ordered: list[ChromosomeRef] = []
        for result in component_results:
            ordered.extend(result.orders.get(species, ()))
        final_order[species] = tuple(ordered)

    final_orientation: dict[ChromosomeRef, int] = {}
    for result in component_results:
        final_orientation.update(result.orientation)

    optimized = LayoutState(
        chromosome_order=final_order,
        chromosome_orientation=final_orientation,
    )
    optimized_score = score_crossings(fixture, optimized)
    component_total = sum(
        result.optimum for result in component_results
    )
    if optimized_score.crossings != component_total:
        raise AssertionError(
            "Residual component optima do not add to the full crossing score: "
            f"components={component_total}, full={optimized_score.crossings}"
        )

    diagnostics = tuple(
        ResidualSolveDiagnostics(
            **{
                **result.diagnostics.__dict__,
                "component_workers": worker_count,
            }
        )
        for result in component_results
    )
    states_evaluated = sum(
        item.table_entries_evaluated for item in diagnostics
    )

    layout = ExactLayoutResult(
        initial_state=initial,
        normalized_state=normalized,
        optimized_state=optimized,
        initial_score=initial_score,
        normalized_score=normalized_score,
        optimized_score=optimized_score,
        states_evaluated=states_evaluated,
        optimality_status="proven optimum",
    )

    return ResidualExactResult(
        layout=layout,
        diagnostics=diagnostics,
        component_workers=worker_count,
    )
