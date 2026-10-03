from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import (
    build_layout_audit,
    exact_optimize_layer_dp,
    load_fixture,
    render_validation_report_html,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class ValidationReportTests(unittest.TestCase):
    def test_report_contains_layout_graph_and_audit_sections(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        result = exact_optimize_layer_dp(fixture)
        audit = build_layout_audit(fixture, result.layout)

        html = render_validation_report_html(
            fixture,
            result.layout,
            solver="exact-layer-dynamic-programming",
            solver_details={
                "component_diagnostics": [
                    item.to_dict() for item in result.diagnostics
                ]
            },
        )

        self.assertIn("SynTangle validation report", html)
        self.assertIn("Tangled input → optimized layout", html)
        self.assertIn("Raw incidence structure", html)
        self.assertIn("Structural projection and hard core", html)
        self.assertIn("proven optimum", html)
        self.assertIn(audit.input_fingerprint, html)
        self.assertIn("<svg", html)

    def test_report_marks_intrinsic_crossing_case(self) -> None:
        fixture = load_fixture(
            FIXTURES / "forbidden_subchromosomal_pretty_solution.json"
        )
        result = exact_optimize_layer_dp(fixture)
        html = render_validation_report_html(
            fixture,
            result.layout,
            solver="exact-layer-dynamic-programming",
        )
        self.assertIn("<th>Optimized crossings</th><td>1</td>", html)


if __name__ == "__main__":
    unittest.main()
