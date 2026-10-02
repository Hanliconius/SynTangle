from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from .model import ChromosomeRef, Fixture


@dataclass(frozen=True)
class IncidenceEdge:
    occurrence_id: str
    chromosome: ChromosomeRef
    homology_id: str

    @property
    def chromosome_node(self) -> str:
        return chromosome_node_id(self.chromosome)

    @property
    def homology_node(self) -> str:
        return homology_node_id(self.homology_id)


@dataclass(frozen=True)
class ComponentSummary:
    node_ids: frozenset[str]
    chromosome_refs: tuple[ChromosomeRef, ...]
    homology_ids: tuple[str, ...]
    incidence_count: int
    cycle_rank: int


@dataclass(frozen=True)
class FixtureAnalysis:
    fixture_id: str
    chromosome_component_count: int
    chromosome_component_sizes: tuple[int, ...]
    incidence_cycle_rank_total: int
    hard_kernel_count: int
    incidence_vertex_count: int
    incidence_edge_count: int
    bridge_count: int
    articulation_point_count: int
    incidence_core_vertex_count: int
    incidence_core_edge_count: int
    biconnected_block_count: int
    cycle_basis_size: int
    components: tuple[ComponentSummary, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "fixture_id": self.fixture_id,
            "chromosome_component_count": self.chromosome_component_count,
            "chromosome_component_sizes": list(self.chromosome_component_sizes),
            "incidence_cycle_rank_total": self.incidence_cycle_rank_total,
            "hard_kernel_count": self.hard_kernel_count,
            "incidence_vertex_count": self.incidence_vertex_count,
            "incidence_edge_count": self.incidence_edge_count,
            "bridge_count": self.bridge_count,
            "articulation_point_count": self.articulation_point_count,
            "incidence_core_vertex_count": self.incidence_core_vertex_count,
            "incidence_core_edge_count": self.incidence_core_edge_count,
            "biconnected_block_count": self.biconnected_block_count,
            "cycle_basis_size": self.cycle_basis_size,
            "components": [
                {
                    "chromosomes": [ref.label for ref in component.chromosome_refs],
                    "homology_ids": list(component.homology_ids),
                    "incidence_count": component.incidence_count,
                    "cycle_rank": component.cycle_rank,
                }
                for component in self.components
            ],
        }


class IncidenceGraph:
    """Sparse chromosome↔homology multigraph.

    Each block occurrence is retained as a distinct incidence edge. Parallel
    occurrences therefore contribute separately to E in cycle-rank accounting,
    while adjacency sets are sufficient for connectivity traversal.
    """

    def __init__(
        self,
        chromosome_nodes: Iterable[ChromosomeRef],
        homology_ids: Iterable[str],
        edges: Iterable[IncidenceEdge],
    ) -> None:
        self.chromosome_refs = tuple(sorted(set(chromosome_nodes)))
        self.homology_ids = tuple(sorted(set(homology_ids)))
        self.edges = tuple(edges)

        self._chromosome_by_node = {
            chromosome_node_id(ref): ref for ref in self.chromosome_refs
        }
        self._homology_by_node = {
            homology_node_id(homology_id): homology_id for homology_id in self.homology_ids
        }
        self._adjacency: dict[str, set[str]] = {
            **{node: set() for node in self._chromosome_by_node},
            **{node: set() for node in self._homology_by_node},
        }

        for edge in self.edges:
            cnode = edge.chromosome_node
            hnode = edge.homology_node
            if cnode not in self._adjacency:
                raise ValueError(f"Incidence edge references unknown chromosome node: {cnode}")
            if hnode not in self._adjacency:
                raise ValueError(f"Incidence edge references unknown homology node: {hnode}")
            self._adjacency[cnode].add(hnode)
            self._adjacency[hnode].add(cnode)

    @property
    def vertex_count(self) -> int:
        return len(self._adjacency)

    @property
    def edge_count(self) -> int:
        return len(self.edges)

    @property
    def node_ids(self) -> tuple[str, ...]:
        return tuple(sorted(self._adjacency))

    @property
    def edge_endpoints(self) -> tuple[tuple[str, str], ...]:
        return tuple(
            (edge.chromosome_node, edge.homology_node)
            for edge in self.edges
        )

    def connected_components(self) -> tuple[frozenset[str], ...]:
        unseen = set(self._adjacency)
        components: list[frozenset[str]] = []

        while unseen:
            seed = min(unseen)
            stack = [seed]
            seen: set[str] = set()
            while stack:
                node = stack.pop()
                if node in seen:
                    continue
                seen.add(node)
                unseen.discard(node)
                stack.extend(sorted(self._adjacency[node] - seen, reverse=True))
            components.append(frozenset(seen))

        components.sort(key=lambda comp: min(comp))
        return tuple(components)

    def summarize_component(self, component: frozenset[str]) -> ComponentSummary:
        chromosome_refs = tuple(
            sorted(self._chromosome_by_node[node] for node in component if node in self._chromosome_by_node)
        )
        homology_ids = tuple(
            sorted(self._homology_by_node[node] for node in component if node in self._homology_by_node)
        )
        incidence_count = sum(
            1
            for edge in self.edges
            if edge.chromosome_node in component and edge.homology_node in component
        )
        # Every value returned by connected_components is connected, so c = 1.
        cycle_rank = incidence_count - len(component) + 1
        if cycle_rank < 0:
            raise AssertionError("Connected incidence component has negative cycle rank")
        return ComponentSummary(
            node_ids=component,
            chromosome_refs=chromosome_refs,
            homology_ids=homology_ids,
            incidence_count=incidence_count,
            cycle_rank=cycle_rank,
        )


def chromosome_node_id(ref: ChromosomeRef) -> str:
    return f"chrom::{ref.species_id}::{ref.chromosome_id}"


def homology_node_id(homology_id: str) -> str:
    return f"homology::{homology_id}"


def build_incidence_graph(fixture: Fixture) -> IncidenceGraph:
    homology_ids: set[str] = set()
    edges: list[IncidenceEdge] = []

    for chromosome in fixture.chromosomes:
        for block in chromosome.blocks:
            homology_ids.add(block.homology_id)
            edges.append(
                IncidenceEdge(
                    occurrence_id=block.occurrence_id,
                    chromosome=chromosome.ref,
                    homology_id=block.homology_id,
                )
            )

    return IncidenceGraph(
        chromosome_nodes=fixture.chromosome_refs,
        homology_ids=homology_ids,
        edges=edges,
    )


def analyze_fixture(fixture: Fixture) -> FixtureAnalysis:
    # Imported here to avoid a module-level circular import:
    # decomposition depends on IncidenceGraph.
    from .decomposition import decompose_incidence_graph

    graph = build_incidence_graph(fixture)
    component_summaries = tuple(
        graph.summarize_component(component) for component in graph.connected_components()
    )

    chromosome_components = tuple(
        component for component in component_summaries if component.chromosome_refs
    )
    chromosome_component_sizes = tuple(
        sorted(len(component.chromosome_refs) for component in chromosome_components)
    )
    cycle_rank_total = sum(component.cycle_rank for component in component_summaries)
    decomposition = decompose_incidence_graph(graph)

    return FixtureAnalysis(
        fixture_id=fixture.fixture_id,
        chromosome_component_count=len(chromosome_components),
        chromosome_component_sizes=chromosome_component_sizes,
        incidence_cycle_rank_total=cycle_rank_total,
        hard_kernel_count=len(decomposition.hard_kernels),
        incidence_vertex_count=graph.vertex_count,
        incidence_edge_count=graph.edge_count,
        bridge_count=len(decomposition.bridge_edge_ids),
        articulation_point_count=len(decomposition.articulation_points),
        incidence_core_vertex_count=len(decomposition.core_node_ids),
        incidence_core_edge_count=len(decomposition.core_edge_ids),
        biconnected_block_count=len(decomposition.biconnected_blocks),
        cycle_basis_size=len(decomposition.cycle_basis),
        components=component_summaries,
    )
