from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import (
    exact_optimize_layer_dp,
    exact_optimize_small,
    load_fixture,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class LayerDynamicProgrammingTests(unittest.TestCase):
    def test_dp_matches_exhaustive_optimum_on_small_fixtures(self) -> None:
        names = (
            "layout_noise_only.json",
            "forbidden_subchromosomal_pretty_solution.json",
            "fusion_chain_tree.json",
            "fusion_chain_closed_cycle.json",
            "balanced_orientation_cycle.json",
        )

        for name in names:
            fixture = load_fixture(FIXTURES / name)
            exhaustive = exact_optimize_small(fixture)
            dp = exact_optimize_layer_dp(fixture)
            self.assertEqual(
                dp.layout.optimized_score.crossings,
                exhaustive.optimized_score.crossings,
                name,
            )
            self.assertEqual(dp.layout.optimality_status, "proven optimum")

    def test_layout_noise_is_removed_by_component_normalization(self) -> None:
        fixture = load_fixture(FIXTURES / "layout_noise_only.json")
        result = exact_optimize_layer_dp(fixture)
        self.assertEqual(result.layout.initial_score.crossings, 3)
        self.assertEqual(result.layout.normalized_score.crossings, 0)
        self.assertEqual(result.layout.optimized_score.crossings, 0)

    def test_orientation_propagation_reduces_assignment_space(self) -> None:
        fixture = load_fixture(FIXTURES / "balanced_orientation_cycle.json")
        result = exact_optimize_layer_dp(fixture)
        self.assertEqual(len(result.diagnostics), 1)
        self.assertEqual(result.diagnostics[0].orientation_assignments, 2)

    def test_perfect_fixture_remains_componentwise_tiny(self) -> None:
        fixture = load_fixture(FIXTURES / "perfect_1to1_30x3.json")
        result = exact_optimize_layer_dp(fixture)
        self.assertEqual(len(result.diagnostics), 30)
        self.assertTrue(
            all(item.chromosome_count == 3 for item in result.diagnostics)
        )
        self.assertEqual(result.layout.optimized_score.crossings, 0)


if __name__ == "__main__":
    unittest.main()
