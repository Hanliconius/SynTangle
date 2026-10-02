from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import build_incidence_graph, fundamental_cycle_basis, load_fixture


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class CycleBasisTests(unittest.TestCase):
    def test_cycle_basis_size_matches_cycle_rank_for_all_fixtures(self) -> None:
        for path in sorted(FIXTURES.glob("*.json")):
            if path.name == "fixture.schema.json":
                continue
            fixture = load_fixture(path)
            graph = build_incidence_graph(fixture)
            basis = fundamental_cycle_basis(graph)
            expected_rank = sum(
                graph.summarize_component(component).cycle_rank
                for component in graph.connected_components()
            )
            self.assertEqual(len(basis), expected_rank, fixture.fixture_id)

    def test_tree_has_empty_basis(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_tree.json")
        graph = build_incidence_graph(fixture)
        self.assertEqual(fundamental_cycle_basis(graph), ())

    def test_closed_cycle_has_one_basis_cycle(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        graph = build_incidence_graph(fixture)
        basis = fundamental_cycle_basis(graph)
        self.assertEqual(len(basis), 1)
        self.assertEqual(len(basis[0].edge_ids), 8)
        self.assertEqual(len(basis[0].node_ids), 8)


if __name__ == "__main__":
    unittest.main()
