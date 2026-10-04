from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    ROOT / "validation" / "benchmark" / "analyze_factorial_benchmark.py"
)

spec = importlib.util.spec_from_file_location(
    "analyze_factorial_benchmark",
    MODULE_PATH,
)
if spec is None or spec.loader is None:
    raise RuntimeError("Could not load factorial benchmark analyzer")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def row(
    case_id: str,
    axis: str,
    *,
    species: int,
    chromosomes: int,
    events: int,
    wall: float,
) -> dict[str, str]:
    return {
        "case_id": case_id,
        "lineage_model": "independent",
        "benchmark_axis": axis,
        "species_count": str(species),
        "ancestor_chromosomes": str(chromosomes),
        "events_per_branch": str(events),
        "wall_seconds": str(wall),
        "optimality_status": "proven optimum",
        "solver": "exact-layer-dp",
        "max_residual_component_variables": "3",
        "residual_treewidth_upper_bound": "2",
        "states_or_nodes_evaluated": "10",
        "optimized_crossings": "0",
    }


class FactorialBenchmarkAnalysisTests(unittest.TestCase):
    def test_axes_are_kept_separate_even_at_shared_coordinates(self) -> None:
        rows = [
            row(
                "species_baseline",
                "species",
                species=8,
                chromosomes=16,
                events=1,
                wall=2.0,
            ),
            row(
                "chromosome_baseline",
                "chromosomes",
                species=8,
                chromosomes=16,
                events=1,
                wall=3.0,
            ),
            row(
                "event_baseline",
                "rearrangements",
                species=8,
                chromosomes=16,
                events=1,
                wall=4.0,
            ),
        ]

        self.assertEqual(
            [item["case_id"] for item in module.axis_rows(rows, "species")],
            ["species_baseline"],
        )
        self.assertEqual(
            [
                item["case_id"]
                for item in module.axis_rows(rows, "chromosomes")
            ],
            ["chromosome_baseline"],
        )
        self.assertEqual(
            [
                item["case_id"]
                for item in module.axis_rows(rows, "rearrangements")
            ],
            ["event_baseline"],
        )

    def test_species_axis_sorts_by_species_count(self) -> None:
        rows = [
            row("s20", "species", species=20, chromosomes=16, events=1, wall=8),
            row("s4", "species", species=4, chromosomes=16, events=1, wall=1),
            row("s12", "species", species=12, chromosomes=16, events=1, wall=4),
        ]
        self.assertEqual(
            [item["case_id"] for item in module.axis_rows(rows, "species")],
            ["s4", "s12", "s20"],
        )

    def test_partial_results_produce_a_valid_section(self) -> None:
        rows = [
            row("e0", "rearrangements", species=8, chromosomes=16, events=0, wall=1),
            row("e2", "rearrangements", species=8, chromosomes=16, events=2, wall=5),
        ]
        rendered = "\n".join(
            module.section(rows, "rearrangements", "Rearrangement burden")
        )
        self.assertIn("e", rendered.lower())
        self.assertIn("5.00x", rendered)


if __name__ == "__main__":
    unittest.main()
