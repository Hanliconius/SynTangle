from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import (
    exact_optimize_layer_dp,
    load_fixture,
    optimize_branch_and_bound,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class BranchAndBoundTests(unittest.TestCase):
    def test_branch_and_bound_matches_exact_dp_on_small_fixtures(self) -> None:
        names = (
            "layout_noise_only.json",
            "forbidden_subchromosomal_pretty_solution.json",
            "fusion_chain_tree.json",
            "fusion_chain_closed_cycle.json",
            "balanced_orientation_cycle.json",
        )

        for name in names:
            fixture = load_fixture(FIXTURES / name)
            exact = exact_optimize_layer_dp(fixture)
            bounded = optimize_branch_and_bound(
                fixture,
                node_cap_per_component=100000,
                local_restarts=3,
                seed=41,
            )
            self.assertEqual(
                bounded.layout.optimized_score.crossings,
                exact.layout.optimized_score.crossings,
                name,
            )
            self.assertEqual(
                bounded.lower_bound,
                bounded.upper_bound,
                name,
            )
            self.assertEqual(
                bounded.layout.optimality_status,
                "proven optimum",
                name,
            )

    def test_cutoff_reports_valid_bound_instead_of_false_optimum(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        exact = exact_optimize_layer_dp(fixture)
        bounded = optimize_branch_and_bound(
            fixture,
            node_cap_per_component=1,
            local_restarts=2,
            seed=43,
        )

        self.assertLessEqual(
            bounded.lower_bound,
            exact.layout.optimized_score.crossings,
        )
        self.assertGreaterEqual(
            bounded.upper_bound,
            exact.layout.optimized_score.crossings,
        )
        self.assertEqual(
            bounded.gap,
            bounded.upper_bound - bounded.lower_bound,
        )
        if bounded.gap > 0:
            self.assertEqual(
                bounded.layout.optimality_status,
                "bounded best known",
            )

    def test_zero_crossing_component_is_immediately_proven(self) -> None:
        fixture = load_fixture(FIXTURES / "layout_noise_only.json")
        bounded = optimize_branch_and_bound(
            fixture,
            node_cap_per_component=1,
            local_restarts=2,
            seed=47,
        )
        self.assertEqual(bounded.upper_bound, 0)
        self.assertEqual(bounded.lower_bound, 0)
        self.assertEqual(bounded.gap, 0)
        self.assertEqual(
            bounded.layout.optimality_status,
            "proven optimum",
        )


if __name__ == "__main__":
    unittest.main()
