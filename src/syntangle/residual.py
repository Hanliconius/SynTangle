from __future__ import annotations

from dataclasses import dataclass

from .incidence import build_incidence_graph, chromosome_node_id
from .model import ChromosomeRef, Fixture
from .orientation_space import orientation_basis


@dataclass(frozen=True)
class ResidualVariable:
    variable_id: str
    kind: str
    component_id: int
    species_id: str | None
    chromosome_refs: tuple[ChromosomeRef, ...]


@dataclass(frozen=True)
class ResidualFactor:
    factor_id: str
    component_id: int
    species_left: str
    species_right: str
    variable_ids: tuple[str, ...]
    homology_ids: tuple[str, ...]


@dataclass(frozen=True)
class ResidualComponentSummary:
    component_id: int
    variable_ids: tuple[str, ...]
    factor_ids: tuple[str, ...]

    @property
    def variable_count(self) -> int:
        return len(self.variable_ids)


@dataclass(frozen=True)
class ResidualFactorization:
    """Variable/factor view of unresolved legal layout decisions.

    This is deliberately distinct from the biological incidence graph.
    Variables are legal display decisions (whole-component chromosome orders
    and GF(2) free orientation groups). Factors are adjacent-species crossing
    terms. The dependency graph is conservative: a factor may list a variable
    that proves irrelevant in a particular state, but it must not omit a
    variable capable of changing that factor. Therefore, if this graph
    disconnects, the current crossing objective exactly factorizes over those
    disconnected residual pieces.

    Stage 18 introduced this factorization and separator diagnostics. Stage 20
    also uses it as an active exact solver representation: evaluated factor
    tables may contract conservative scopes, and the remaining graph is
    recursively split/eliminated by the residual factor solver.
    """

    variables: tuple[ResidualVariable, ...]
    factors: tuple[ResidualFactor, ...]
    objective_components: tuple[ResidualComponentSummary, ...]
    isolated_variable_ids: tuple[str, ...]
    articulation_variable_ids: tuple[str, ...]
    min_fill_treewidth_upper_bound: int

    @property
    def objective_variable_count(self) -> int:
        return sum(item.variable_count for item in self.objective_components)

    @property
    def objective_component_count(self) -> int:
        return len(self.objective_components)

    @property
    def max_objective_component_variables(self) -> int:
        return max(
            (item.variable_count for item in self.objective_components),
            default=0,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "variable_count": len(self.variables),
            "factor_count": len(self.factors),
            "objective_variable_count": self.objective_variable_count,
            "objective_component_count": self.objective_component_count,
            "max_objective_component_variables": (
                self.max_objective_component_variables
            ),
            "isolated_variable_count": len(self.isolated_variable_ids),
            "articulation_variable_count": len(
                self.articulation_variable_ids
            ),
            "min_fill_treewidth_upper_bound": (
                self.min_fill_treewidth_upper_bound
            ),
            "isolated_variable_ids": list(self.isolated_variable_ids),
            "articulation_variable_ids": list(
                self.articulation_variable_ids
            ),
            "objective_components": [
                {
                    "component_id": item.component_id,
                    "variable_ids": list(item.variable_ids),
                    "factor_ids": list(item.factor_ids),
                }
                for item in self.objective_components
            ],
        }


def _primal_adjacency(
    variable_ids: tuple[str, ...],
    factors: tuple[ResidualFactor, ...],
) -> dict[str, set[str]]:
    adjacency = {variable_id: set() for variable_id in variable_ids}
    for factor in factors:
        ids = factor.variable_ids
        for index, a in enumerate(ids):
            for b in ids[index + 1 :]:
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


def _min_fill_treewidth_upper_bound(
    adjacency: dict[str, set[str]],
) -> int:
    """Greedy min-fill elimination width.

    This is an upper bound from one elimination heuristic, not an exact
    treewidth calculation. It is recorded as a scaling diagnostic only.
    """

    graph = {
        node: set(neighbors)
        for node, neighbors in adjacency.items()
    }
    width = 0

    while graph:
        def key(node: str) -> tuple[int, int, str]:
            neighbors = sorted(graph[node])
            missing = 0
            for index, a in enumerate(neighbors):
                for b in neighbors[index + 1 :]:
                    if b not in graph[a]:
                        missing += 1
            return (missing, len(neighbors), node)

        node = min(graph, key=key)
        neighbors = list(graph[node])
        width = max(width, len(neighbors))

        for index, a in enumerate(neighbors):
            for b in neighbors[index + 1 :]:
                graph[a].add(b)
                graph[b].add(a)

        for neighbor in neighbors:
            graph[neighbor].discard(node)
        del graph[node]

    return width


def build_residual_factorization(
    fixture: Fixture,
) -> ResidualFactorization:
    """Build the residual decision-factor graph for the current objective."""

    graph = build_incidence_graph(fixture)
    incidence_components = graph.connected_components()

    component_of: dict[ChromosomeRef, int] = {}
    for component_id, nodes in enumerate(incidence_components):
        for ref in fixture.chromosome_refs:
            if chromosome_node_id(ref) in nodes:
                component_of[ref] = component_id

    homology_refs: dict[str, set[ChromosomeRef]] = {}
    for chromosome in fixture.chromosomes:
        for block in chromosome.blocks:
            homology_refs.setdefault(block.homology_id, set()).add(
                chromosome.ref
            )

    variables: list[ResidualVariable] = []
    factors: list[ResidualFactor] = []

    for component_id, component_nodes in enumerate(incidence_components):
        refs = tuple(
            sorted(
                ref
                for ref in fixture.chromosome_refs
                if component_of[ref] == component_id
            )
        )
        ref_set = set(refs)

        refs_by_species = {
            species: tuple(
                ref for ref in refs if ref.species_id == species
            )
            for species in fixture.species_ids
        }

        order_variable_by_species: dict[str, str] = {}
        for species, species_refs in refs_by_species.items():
            if len(species_refs) <= 1:
                continue
            variable_id = f"order::{component_id}::{species}"
            order_variable_by_species[species] = variable_id
            variables.append(
                ResidualVariable(
                    variable_id=variable_id,
                    kind="order",
                    component_id=component_id,
                    species_id=species,
                    chromosome_refs=species_refs,
                )
            )

        basis = orientation_basis(fixture, refs)
        orientation_group_by_ref: dict[ChromosomeRef, str] = {}
        for group_index, group in enumerate(basis.free_flip_groups):
            variable_id = f"orient::{component_id}::{group_index}"
            variables.append(
                ResidualVariable(
                    variable_id=variable_id,
                    kind="orientation",
                    component_id=component_id,
                    species_id=None,
                    chromosome_refs=tuple(group),
                )
            )
            for ref in group:
                orientation_group_by_ref[ref] = variable_id

        component_homology_ids = tuple(
            sorted(
                homology_id
                for homology_id, hom_refs in homology_refs.items()
                if hom_refs & ref_set
            )
        )

        for species_left, species_right in zip(
            fixture.species_ids,
            fixture.species_ids[1:],
        ):
            shared_homology: list[str] = []
            participating_refs: set[ChromosomeRef] = set()

            for homology_id in component_homology_ids:
                hom_refs = homology_refs[homology_id] & ref_set
                left_refs = {
                    ref
                    for ref in hom_refs
                    if ref.species_id == species_left
                }
                right_refs = {
                    ref
                    for ref in hom_refs
                    if ref.species_id == species_right
                }
                if not left_refs or not right_refs:
                    continue
                shared_homology.append(homology_id)
                participating_refs.update(left_refs)
                participating_refs.update(right_refs)

            if not shared_homology:
                continue

            variable_ids: set[str] = set()
            left_order = order_variable_by_species.get(species_left)
            right_order = order_variable_by_species.get(species_right)
            if left_order is not None:
                variable_ids.add(left_order)
            if right_order is not None:
                variable_ids.add(right_order)

            for ref in participating_refs:
                variable_ids.add(orientation_group_by_ref[ref])

            factors.append(
                ResidualFactor(
                    factor_id=(
                        f"cross::{component_id}::"
                        f"{species_left}::{species_right}"
                    ),
                    component_id=component_id,
                    species_left=species_left,
                    species_right=species_right,
                    variable_ids=tuple(sorted(variable_ids)),
                    homology_ids=tuple(sorted(shared_homology)),
                )
            )

    variable_ids = tuple(
        variable.variable_id for variable in variables
    )
    active_ids = {
        variable_id
        for factor in factors
        for variable_id in factor.variable_ids
    }
    isolated = tuple(sorted(set(variable_ids) - active_ids))

    factor_by_variable: dict[str, set[str]] = {
        variable_id: set() for variable_id in active_ids
    }
    variables_by_factor = {
        factor.factor_id: set(factor.variable_ids)
        for factor in factors
        if factor.variable_ids
    }
    for factor_id, ids in variables_by_factor.items():
        for variable_id in ids:
            factor_by_variable[variable_id].add(factor_id)

    unseen = set(active_ids)
    objective_components: list[ResidualComponentSummary] = []
    residual_component_id = 0

    while unseen:
        seed = min(unseen)
        stack_variables = [seed]
        seen_variables: set[str] = set()
        seen_factors: set[str] = set()

        while stack_variables:
            variable_id = stack_variables.pop()
            if variable_id in seen_variables:
                continue
            seen_variables.add(variable_id)
            unseen.discard(variable_id)

            for factor_id in factor_by_variable[variable_id]:
                if factor_id in seen_factors:
                    continue
                seen_factors.add(factor_id)
                for neighbor in variables_by_factor[factor_id]:
                    if neighbor not in seen_variables:
                        stack_variables.append(neighbor)

        objective_components.append(
            ResidualComponentSummary(
                component_id=residual_component_id,
                variable_ids=tuple(sorted(seen_variables)),
                factor_ids=tuple(sorted(seen_factors)),
            )
        )
        residual_component_id += 1

    objective_components.sort(
        key=lambda item: (
            -item.variable_count,
            item.variable_ids,
        )
    )
    objective_components = [
        ResidualComponentSummary(
            component_id=index,
            variable_ids=item.variable_ids,
            factor_ids=item.factor_ids,
        )
        for index, item in enumerate(objective_components)
    ]

    active_variable_ids = tuple(sorted(active_ids))
    active_factors = tuple(
        factor
        for factor in factors
        if factor.variable_ids
    )
    primal = _primal_adjacency(active_variable_ids, active_factors)

    return ResidualFactorization(
        variables=tuple(variables),
        factors=tuple(factors),
        objective_components=tuple(objective_components),
        isolated_variable_ids=isolated,
        articulation_variable_ids=_articulation_points(primal),
        min_fill_treewidth_upper_bound=(
            _min_fill_treewidth_upper_bound(primal)
        ),
    )
