from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import (
    exact_optimize_small,
    initial_layout_state,
    load_fixture,
    score_crossings,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class ExactLayoutTests(unittest.TestCase):
    def test_layout_noise_fixture_recovers_expected_zero(self) -> None:
        fixture = load_fixture(FIXTURES / "layout_noise_only.json")
        initial = initial_layout_state(fixture)
        self.assertEqual(
            score_crossings(fixture, initial).crossings,
            fixture.expected["initial_whole_chromosome_crossings"],
        )

        result = exact_optimize_small(fixture)
        self.assertEqual(
            result.optimized_score.crossings,
            fixture.expected["constrained_minimum_crossings"],
        )
        self.assertEqual(
            result.excess_layout_crossings_initial,
            fixture.expected["excess_layout_crossings_initial"],
        )
        self.assertEqual(result.normalized_score.crossings, 0)
        self.assertEqual(result.optimality_status, "proven optimum")

    def test_forbidden_internal_swap_retains_intrinsic_crossing(self) -> None:
        fixture = load_fixture(
            FIXTURES / "forbidden_subchromosomal_pretty_solution.json"
        )
        result = exact_optimize_small(fixture)

        # Whole-chromosome flips are legal, but no internal B/C swap is.
        self.assertEqual(result.optimized_score.crossings, 1)
        self.assertGreater(result.optimized_score.crossings, 0)

        original_orders = {
            chrom.ref.label: tuple(block.homology_id for block in chrom.blocks)
            for chrom in fixture.chromosomes
        }
        self.assertEqual(original_orders["sp1:X"], ("HA", "HB", "HC", "HD"))
        self.assertEqual(original_orders["sp2:Y"], ("HA", "HC", "HB", "HD"))

    def test_balanced_orientation_constraints_are_respected(self) -> None:
        fixture = load_fixture(FIXTURES / "balanced_orientation_cycle.json")
        result = exact_optimize_small(fixture)
        orientation = result.optimized_state.chromosome_orientation

        for constraint in fixture.orientation_constraints:
            abit = 0 if orientation[constraint.a] == 1 else 1
            bbit = 0 if orientation[constraint.b] == 1 else 1
            self.assertEqual(abit ^ bbit, constraint.xor)

    def test_solver_changes_only_whole_chromosome_state(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        before = {
            chrom.ref: tuple(
                (block.occurrence_id, block.start, block.end, block.homology_id)
                for block in chrom.blocks
            )
            for chrom in fixture.chromosomes
        }
        result = exact_optimize_small(fixture)
        after = {
            chrom.ref: tuple(
                (block.occurrence_id, block.start, block.end, block.homology_id)
                for block in chrom.blocks
            )
            for chrom in fixture.chromosomes
        }
        self.assertEqual(before, after)
        self.assertEqual(
            set(result.optimized_state.chromosome_orientation),
            set(fixture.chromosome_refs),
        )


if __name__ == "__main__":
    unittest.main()
