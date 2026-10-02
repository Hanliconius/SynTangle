from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import (
    build_layout_audit,
    exact_optimize_small,
    fixture_fingerprint,
    load_fixture,
    render_layout_comparison_svg,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class AuditAndVisualizationTests(unittest.TestCase):
    def test_layout_noise_audit_records_improvement(self) -> None:
        fixture = load_fixture(FIXTURES / "layout_noise_only.json")
        result = exact_optimize_small(fixture)
        audit = build_layout_audit(fixture, result)

        self.assertEqual(audit.objective_before, 3)
        self.assertEqual(audit.objective_after, 0)
        self.assertEqual(audit.excess_layout_crossings_initial, 3)
        self.assertEqual(audit.optimality_status, "proven optimum")
        self.assertGreater(len(audit.whole_chromosome_moves), 0)
        self.assertEqual(audit.input_fingerprint, fixture_fingerprint(fixture))
        self.assertEqual(len(audit.input_fingerprint), 64)

    def test_forbidden_fixture_audit_retains_intrinsic_crossing(self) -> None:
        fixture = load_fixture(
            FIXTURES / "forbidden_subchromosomal_pretty_solution.json"
        )
        result = exact_optimize_small(fixture)
        audit = build_layout_audit(fixture, result)
        self.assertEqual(audit.objective_after, 1)

    def test_svg_contains_both_states_and_crossing_counts(self) -> None:
        fixture = load_fixture(FIXTURES / "layout_noise_only.json")
        result = exact_optimize_small(fixture)
        svg = render_layout_comparison_svg(fixture, result)

        self.assertIn("<svg", svg)
        self.assertIn("Initial — crossings: 3", svg)
        self.assertIn("Optimized — crossings: 0", svg)
        self.assertIn("sp1", svg)
        self.assertIn("sp2", svg)


if __name__ == "__main__":
    unittest.main()
