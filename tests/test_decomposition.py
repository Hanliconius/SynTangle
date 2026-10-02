from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import build_incidence_graph, decompose_incidence_graph, load_fixture


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class GraphDecompositionTests(unittest.TestCase):
    def test_tree_has_only_bridges_and_empty_two_core(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_tree.json")
        graph = build_incidence_graph(fixture)
        result = decompose_incidence_graph(graph)
        self.assertEqual(len(result.bridge_edge_ids), graph.edge_count)
        self.assertEqual(len(result.core_node_ids), 0)
        self.assertEqual(len(result.core_edge_ids), 0)
        self.assertEqual(len(result.hard_kernels), 0)
        self.assertEqual(len(result.biconnected_blocks), graph.edge_count)

    def test_closed_cycle_reduces_to_one_eight_node_kernel(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        graph = build_incidence_graph(fixture)
        result = decompose_incidence_graph(graph)
        self.assertEqual(len(result.bridge_edge_ids), 4)
        self.assertEqual(len(result.articulation_points), 4)
        self.assertEqual(len(result.core_node_ids), 8)
        self.assertEqual(len(result.core_edge_ids), 8)
        self.assertEqual(len(result.hard_kernels), 1)
        self.assertEqual(len(result.hard_kernels[0]), 8)
        self.assertEqual(len(result.biconnected_blocks), 5)

    def test_perfect_one_to_one_has_no_hard_core(self) -> None:
        fixture = load_fixture(FIXTURES / "perfect_1to1_30x3.json")
        graph = build_incidence_graph(fixture)
        result = decompose_incidence_graph(graph)
        self.assertEqual(len(result.core_node_ids), 0)
        self.assertEqual(len(result.hard_kernels), 0)


if __name__ == "__main__":
    unittest.main()
