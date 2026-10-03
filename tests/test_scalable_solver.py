from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import load_fixture, optimize_auto, optimize_local_search


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class ScalableSolverTests(unittest.TestCase):
    def test_local_search_removes_layout_noise(self) -> None:
        fixture = load_fixture(FIXTURES / "layout_noise_only.json")
        result = optimize_local_search(fixture, restarts=3, seed=19)
        self.assertEqual(result.layout.optimized_score.crossings, 0)
        self.assertEqual(result.layout.optimality_status, "best known")

    def test_local_search_cannot_use_forbidden_internal_swap(self) -> None:
        fixture = load_fixture(
            FIXTURES / "forbidden_subchromosomal_pretty_solution.json"
        )
        result = optimize_local_search(fixture, restarts=4, seed=23)
        self.assertEqual(result.layout.optimized_score.crossings, 1)

    def test_local_search_preserves_hard_orientation_equations(self) -> None:
        fixture = load_fixture(FIXTURES / "balanced_orientation_cycle.json")
        result = optimize_local_search(fixture, restarts=4, seed=29)
        orientation = result.layout.optimized_state.chromosome_orientation
        for constraint in fixture.orientation_constraints:
            abit = 0 if orientation[constraint.a] == 1 else 1
            bbit = 0 if orientation[constraint.b] == 1 else 1
            self.assertEqual(abit ^ bbit, constraint.xor)

    def test_auto_prefers_exact_solver_when_feasible(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        result = optimize_auto(fixture)
        self.assertEqual(result.solver, "exact-layer-dynamic-programming")
        self.assertEqual(result.layout.optimality_status, "proven optimum")

    def test_auto_routes_oversized_exact_case_to_branch_and_bound(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        result = optimize_auto(
            fixture,
            transition_cap_per_component=1,
            branch_node_cap_per_component=100000,
            local_restarts=3,
            seed=31,
        )
        self.assertEqual(result.solver, "monotone-component-branch-and-bound")
        self.assertIn(
            result.layout.optimality_status,
            {"proven optimum", "bounded best known"},
        )
        self.assertIn("fallback_reason", result.details)
        self.assertIn("lower_bound", result.details)
        self.assertIn("upper_bound", result.details)

    def test_seeded_local_search_is_reproducible(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        first = optimize_local_search(fixture, restarts=5, seed=37)
        second = optimize_local_search(fixture, restarts=5, seed=37)
        self.assertEqual(
            first.layout.optimized_state.to_dict(),
            second.layout.optimized_state.to_dict(),
        )
        self.assertEqual(
            first.layout.optimized_score.crossings,
            second.layout.optimized_score.crossings,
        )


if __name__ == "__main__":
    unittest.main()
