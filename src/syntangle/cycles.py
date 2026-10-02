from __future__ import annotations

from dataclasses import dataclass
from collections import deque

from .incidence import IncidenceGraph


@dataclass(frozen=True)
class CycleBasisElement:
    edge_ids: tuple[int, ...]
    node_ids: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "edge_ids": list(self.edge_ids),
            "node_ids": list(self.node_ids),
        }


class _DisjointSet:
    def __init__(self, nodes: tuple[str, ...]) -> None:
        self.parent = {node: node for node in nodes}
        self.rank = {node: 0 for node in nodes}

    def find(self, node: str) -> str:
        parent = self.parent[node]
        if parent != node:
            self.parent[node] = self.find(parent)
        return self.parent[node]

    def union(self, a: str, b: str) -> bool:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True


def _tree_path_edges(
    start: str,
    goal: str,
    tree_adjacency: dict[str, list[tuple[str, int]]],
) -> tuple[int, ...]:
    queue: deque[str] = deque([start])
    parent: dict[str, tuple[str, int] | None] = {start: None}

    while queue:
        node = queue.popleft()
        if node == goal:
            break
        for neighbor, edge_id in tree_adjacency[node]:
            if neighbor in parent:
                continue
            parent[neighbor] = (node, edge_id)
            queue.append(neighbor)

    if goal not in parent:
        raise AssertionError("Spanning forest path missing for non-tree edge")

    path: list[int] = []
    node = goal
    while node != start:
        prev = parent[node]
        assert prev is not None
        parent_node, edge_id = prev
        path.append(edge_id)
        node = parent_node
    path.reverse()
    return tuple(path)


def fundamental_cycle_basis(
    graph: IncidenceGraph,
    *,
    node_ids: frozenset[str] | None = None,
    edge_ids: tuple[int, ...] | None = None,
) -> tuple[CycleBasisElement, ...]:
    """Return a deterministic fundamental cycle basis.

    The algorithm builds a spanning forest using edge IDs, then adds one basis
    cycle for every non-tree edge. It therefore runs in polynomial time and
    never enumerates all simple cycles.

    Parallel incidence edges are handled correctly: one can enter the spanning
    forest and another then closes a two-edge multigraph cycle.
    """

    selected_nodes = set(graph.node_ids if node_ids is None else node_ids)
    selected_edges = tuple(
        range(graph.edge_count) if edge_ids is None else edge_ids
    )

    for edge_id in selected_edges:
        a, b = graph.edge_endpoints[edge_id]
        if a not in selected_nodes or b not in selected_nodes:
            raise ValueError("Selected edge has endpoint outside selected node set")

    dsu = _DisjointSet(tuple(sorted(selected_nodes)))
    tree_adjacency = {node: [] for node in selected_nodes}
    non_tree_edges: list[int] = []

    for edge_id in sorted(selected_edges):
        a, b = graph.edge_endpoints[edge_id]
        if dsu.union(a, b):
            tree_adjacency[a].append((b, edge_id))
            tree_adjacency[b].append((a, edge_id))
        else:
            non_tree_edges.append(edge_id)

    basis: list[CycleBasisElement] = []
    for edge_id in non_tree_edges:
        a, b = graph.edge_endpoints[edge_id]
        tree_path = _tree_path_edges(a, b, tree_adjacency)
        cycle_edges = tuple(sorted((*tree_path, edge_id)))
        cycle_nodes: set[str] = set()
        for cycle_edge_id in cycle_edges:
            u, v = graph.edge_endpoints[cycle_edge_id]
            cycle_nodes.add(u)
            cycle_nodes.add(v)
        basis.append(
            CycleBasisElement(
                edge_ids=cycle_edges,
                node_ids=tuple(sorted(cycle_nodes)),
            )
        )

    basis.sort(key=lambda cycle: (cycle.edge_ids, cycle.node_ids))
    return tuple(basis)
