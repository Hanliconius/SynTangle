from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import (
    build_residual_factorization,
    load_fixture,
    optimize_branch_and_bound,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class ResidualFactorizationTests(unittest.TestCase):
    def test_perfect_1to1_retains_independent_residual_pieces(self) -> None:
        fixture = load_fixture(FIXTURES / "perfect_1to1_30x3.json")
        residual = build_residual_factorization(fixture)

        self.assertGreater(len(residual.variables), 0)
        self.assertGreater(len(residual.factors), 0)
        self.assertGreaterEqual(residual.objective_component_count, 30)
        self.assertLessEqual(
            residual.max_objective_component_variables,
            residual.objective_variable_count,
        )

        # Residual factorization must never reconnect different exact incidence
        # components. Variable IDs carry their parent component ID.
        for component in residual.objective_components:
            parent_ids = {
                variable_id.split("::")[1]
                for variable_id in component.variable_ids
            }
            self.assertEqual(len(parent_ids), 1)

    def test_treewidth_diagnostic_is_a_nonnegative_upper_bound(self) -> None:
        fixture = load_fixture(
            FIXTURES / "fusion_chain_closed_cycle.json"
        )
        residual = build_residual_factorization(fixture)
        self.assertGreaterEqual(
            residual.min_fill_treewidth_upper_bound,
            0,
        )
        self.assertLessEqual(
            residual.min_fill_treewidth_upper_bound,
            max(residual.objective_variable_count - 1, 0),
        )

    def test_parallel_incidence_components_match_serial_result(self) -> None:
        fixture = load_fixture(FIXTURES / "perfect_1to1_30x3.json")
        serial = optimize_branch_and_bound(
            fixture,
            node_cap_per_component=1000,
            local_restarts=1,
            seed=71,
            component_workers=1,
        )
        parallel = optimize_branch_and_bound(
            fixture,
            node_cap_per_component=1000,
            local_restarts=1,
            seed=71,
            component_workers=2,
        )

        self.assertEqual(
            parallel.layout.optimized_score.crossings,
            serial.layout.optimized_score.crossings,
        )
        self.assertEqual(parallel.lower_bound, serial.lower_bound)
        self.assertEqual(parallel.upper_bound, serial.upper_bound)
        self.assertEqual(parallel.gap, serial.gap)
        self.assertEqual(parallel.component_workers, 2)


if __name__ == "__main__":
    unittest.main()
