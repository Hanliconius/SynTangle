from __future__ import annotations

import copy
import json
import unittest
from pathlib import Path

from syntangle import FixtureValidationError, analyze_fixture, fixture_from_dict, load_fixture


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class Stage0FixtureTests(unittest.TestCase):
    def test_expected_incidence_summaries(self) -> None:
        for path in sorted(FIXTURES.glob("*.json")):
            if path.name == "fixture.schema.json":
                continue
            fixture = load_fixture(path)
            analysis = analyze_fixture(fixture)
            expected = fixture.expected

            if "chromosome_component_count" in expected:
                self.assertEqual(
                    analysis.chromosome_component_count,
                    expected["chromosome_component_count"],
                    fixture.fixture_id,
                )
            if "chromosome_component_sizes" in expected:
                self.assertEqual(
                    list(analysis.chromosome_component_sizes),
                    sorted(expected["chromosome_component_sizes"]),
                    fixture.fixture_id,
                )
            if "incidence_cycle_rank_total" in expected:
                self.assertEqual(
                    analysis.incidence_cycle_rank_total,
                    expected["incidence_cycle_rank_total"],
                    fixture.fixture_id,
                )
            if "hard_kernel_count" in expected:
                self.assertEqual(
                    analysis.hard_kernel_count,
                    expected["hard_kernel_count"],
                    fixture.fixture_id,
                )

    def test_first_milestone_exactly(self) -> None:
        expected = {
            "perfect_1to1_30x3.json": (30, 0, 0),
            "fusion_chain_tree.json": (1, 0, 0),
            "fusion_chain_closed_cycle.json": (1, 1, 1),
        }
        for filename, (components, cycle_rank, hard_kernels) in expected.items():
            result = analyze_fixture(load_fixture(FIXTURES / filename))
            self.assertEqual(result.chromosome_component_count, components)
            self.assertEqual(result.incidence_cycle_rank_total, cycle_rank)
            self.assertEqual(result.hard_kernel_count, hard_kernels)

    def test_tree_fixture_incidence_counts(self) -> None:
        result = analyze_fixture(load_fixture(FIXTURES / "fusion_chain_tree.json"))
        self.assertEqual(result.incidence_vertex_count, 11)
        self.assertEqual(result.incidence_edge_count, 10)

    def test_duplicate_occurrence_ids_are_rejected(self) -> None:
        path = FIXTURES / "fusion_chain_tree.json"
        data = json.loads(path.read_text())
        broken = copy.deepcopy(data)
        first = broken["species"][0]["chromosomes"][0]["blocks"][0]["occurrence_id"]
        broken["species"][1]["chromosomes"][0]["blocks"][0]["occurrence_id"] = first
        with self.assertRaises(FixtureValidationError):
            fixture_from_dict(broken)

    def test_block_order_is_coordinate_canonical(self) -> None:
        data = json.loads((FIXTURES / "fusion_chain_tree.json").read_text())
        blocks = data["species"][0]["chromosomes"][0]["blocks"]
        data["species"][0]["chromosomes"][0]["blocks"] = list(reversed(blocks))
        fixture = fixture_from_dict(data)
        chromosome = next(c for c in fixture.chromosomes if c.ref.label == "sp1:A")
        self.assertEqual([b.homology_id for b in chromosome.blocks], ["H1", "H30"])


if __name__ == "__main__":
    unittest.main()
