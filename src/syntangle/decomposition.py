from __future__ import annotations

from dataclasses import dataclass

from .incidence import IncidenceGraph


@dataclass(frozen=True)
class BiconnectedBlock:
    node_ids: frozenset[str]
    edge_ids: tuple[int, ...]

    @property
    def edge_count(self) -> int:
        return len(self.edge_ids)


@dataclass(frozen=True)
class GraphDecomposition:
    bridge_edge_ids: tuple[int, ...]
    articulation_points: tuple[str, ...]
    core_node_ids: frozenset[str]
    core_edge_ids: tuple[int, ...]
    biconnected_blocks: tuple[BiconnectedBlock, ...]
    hard_kernels: tuple[frozenset[str], ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "bridge_edge_ids": list(self.bridge_edge_ids),
            "articulation_points": list(self.articulation_points),
            "core_node_ids": sorted(self.core_node_ids),
            "core_edge_ids": list(self.core_edge_ids),
            "biconnected_blocks": [
                {
                    "node_ids": sorted(block.node_ids),
                    "edge_ids": list(block.edge_ids),
                    "edge_count": block.edge_count,
                }
                for block in self.biconnected_blocks
            ],
            "hard_kernels": [sorted(kernel) for kernel in self.hard_kernels],
        }


def _edge_adjacency(graph: IncidenceGraph) -> dict[str, list[tuple[str, int]]]:
    adjacency = {node: [] for node in graph.node_ids}
    for edge_id, (a, b) in enumerate(graph.edge_endpoints):
        adjacency[a].append((b, edge_id))
        adjacency[b].append((a, edge_id))
    for node in adjacency:
        adjacency[node].sort(key=lambda item: (item[0], item[1]))
    return adjacency


def _tarjan_blocks(
    graph: IncidenceGraph,
) -> tuple[tuple[int, ...], tuple[str, ...], tuple[BiconnectedBlock, ...]]:
    adjacency = _edge_adjacency(graph)
    discovery: dict[str, int] = {}
    low: dict[str, int] = {}
    parent_edge: dict[str, int | None] = {}
    edge_stack: list[int] = []
    bridges: set[int] = set()
    articulation: set[str] = set()
    blocks: list[BiconnectedBlock] = []
    time = 0

    def pop_block(stop_edge: int) -> None:
        edge_ids: list[int] = []
        node_ids: set[str] = set()
        while edge_stack:
            edge_id = edge_stack.pop()
            edge_ids.append(edge_id)
            a, b = graph.edge_endpoints[edge_id]
            node_ids.add(a)
            node_ids.add(b)
            if edge_id == stop_edge:
                break
        if edge_ids:
            blocks.append(
                BiconnectedBlock(
                    node_ids=frozenset(node_ids),
                    edge_ids=tuple(sorted(edge_ids)),
                )
            )

    def dfs(node: str, root: str) -> None:
        nonlocal time
        time += 1
        discovery[node] = low[node] = time
        children = 0

        for neighbor, edge_id in adjacency[node]:
            if edge_id == parent_edge.get(node):
                continue

            if neighbor not in discovery:
                parent_edge[neighbor] = edge_id
                children += 1
                edge_stack.append(edge_id)
                dfs(neighbor, root)
                low[node] = min(low[node], low[neighbor])

                if low[neighbor] > discovery[node]:
                    bridges.add(edge_id)

                if low[neighbor] >= discovery[node]:
                    if node != root or children > 1:
                        articulation.add(node)
                    pop_block(edge_id)

            elif discovery[neighbor] < discovery[node]:
                # Back edge to an ancestor. Distinct edge IDs make this work
                # correctly for parallel edges in the incidence multigraph.
                edge_stack.append(edge_id)
                low[node] = min(low[node], discovery[neighbor])

        if node == root and children <= 1:
            articulation.discard(node)

    for root in graph.node_ids:
        if root in discovery:
            continue
        parent_edge[root] = None
        dfs(root, root)
        if edge_stack:
            # Defensive flush for a component whose final edges were not popped.
            pop_block(edge_stack[0])

    blocks.sort(key=lambda block: (min(block.node_ids), block.edge_ids))
    return tuple(sorted(bridges)), tuple(sorted(articulation)), tuple(blocks)


def _two_core(graph: IncidenceGraph) -> tuple[frozenset[str], tuple[int, ...]]:
    adjacency = _edge_adjacency(graph)
    active_nodes = set(graph.node_ids)
    active_edges = set(range(graph.edge_count))
    degree = {node: len(adjacency[node]) for node in graph.node_ids}
    queue = sorted(node for node, deg in degree.items() if deg < 2)

    while queue:
        node = queue.pop(0)
        if node not in active_nodes:
            continue
        active_nodes.remove(node)
        for neighbor, edge_id in adjacency[node]:
            if edge_id not in active_edges:
                continue
            active_edges.remove(edge_id)
            if neighbor in active_nodes:
                degree[neighbor] -= 1
                if degree[neighbor] < 2:
                    queue.append(neighbor)
        queue.sort()

    return frozenset(active_nodes), tuple(sorted(active_edges))


def _induced_components(
    graph: IncidenceGraph, node_ids: frozenset[str], edge_ids: tuple[int, ...]
) -> tuple[frozenset[str], ...]:
    if not node_ids:
        return ()
    adjacency = {node: set() for node in node_ids}
    for edge_id in edge_ids:
        a, b = graph.edge_endpoints[edge_id]
        if a in node_ids and b in node_ids:
            adjacency[a].add(b)
            adjacency[b].add(a)

    unseen = set(node_ids)
    components: list[frozenset[str]] = []
    while unseen:
        root = min(unseen)
        stack = [root]
        seen: set[str] = set()
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            unseen.discard(node)
            stack.extend(sorted(adjacency[node] - seen, reverse=True))
        components.append(frozenset(seen))
    components.sort(key=lambda comp: min(comp))
    return tuple(components)


def decompose_incidence_graph(graph: IncidenceGraph) -> GraphDecomposition:
    bridges, articulation, blocks = _tarjan_blocks(graph)
    core_nodes, core_edges = _two_core(graph)
    hard_kernels = _induced_components(graph, core_nodes, core_edges)
    return GraphDecomposition(
        bridge_edge_ids=bridges,
        articulation_points=articulation,
        core_node_ids=core_nodes,
        core_edge_ids=core_edges,
        biconnected_blocks=blocks,
        hard_kernels=hard_kernels,
    )
