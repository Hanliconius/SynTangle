from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import load_fixture, solve_orientation_constraints


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class OrientationConstraintTests(unittest.TestCase):
    def test_balanced_cycle_propagates_exactly(self) -> None:
        fixture = load_fixture(FIXTURES / "balanced_orientation_cycle.json")
        result = solve_orientation_constraints(fixture)
        self.assertTrue(result.balanced)
        self.assertTrue(result.exact)
        self.assertEqual(result.frustration_index, fixture.expected["frustration_index"])
        self.assertEqual(result.free_bits, fixture.expected["orientation_free_bits"])
        self.assertEqual(len(result.unsatisfied), 0)

        for constraint in fixture.orientation_constraints:
            observed = result.assignment[constraint.a] ^ result.assignment[constraint.b]
            self.assertEqual(observed, constraint.xor)

    def test_frustrated_cycle_has_exact_frustration_one(self) -> None:
        fixture = load_fixture(FIXTURES / "frustrated_orientation_cycle.json")
        result = solve_orientation_constraints(fixture)
        self.assertFalse(result.balanced)
        self.assertTrue(result.exact)
        self.assertEqual(result.frustration_index, fixture.expected["frustration_index"])
        self.assertEqual(len(result.unsatisfied), 1)

    def test_global_reversal_symmetry_preserves_equations(self) -> None:
        fixture = load_fixture(FIXTURES / "balanced_orientation_cycle.json")
        result = solve_orientation_constraints(fixture)
        flipped = {ref: bit ^ 1 for ref, bit in result.assignment.items()}
        for constraint in fixture.orientation_constraints:
            observed = flipped[constraint.a] ^ flipped[constraint.b]
            self.assertEqual(observed, constraint.xor)

    def test_unconstrained_chromosomes_are_free_bits(self) -> None:
        fixture = load_fixture(FIXTURES / "layout_noise_only.json")
        result = solve_orientation_constraints(fixture)
        self.assertTrue(result.balanced)
        self.assertEqual(result.frustration_index, 0)
        self.assertEqual(result.free_bits, len(fixture.chromosomes))


if __name__ == "__main__":
    unittest.main()
