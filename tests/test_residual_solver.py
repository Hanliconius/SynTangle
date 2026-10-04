from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import (
    exact_optimize_layer_dp,
    exact_optimize_residual_factor_graph,
    fixture_from_dict,
    load_fixture,
    optimize_auto,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


def orientation_chain_fixture(
    n_species: int,
    *,
    n_blocks: int = 3,
):
    species = []
    for index in range(n_species):
        species_id = f"sp{index + 1}"
        blocks = []
        for block_index in range(n_blocks):
            start = block_index * 100
            blocks.append(
                {
                    "occurrence_id": (
                        f"{species_id}:chr1:H{block_index + 1}"
                    ),
                    "homology_id": f"H{block_index + 1}",
                    "start": start,
                    "end": start + 90,
                    "strand": "+",
                }
            )
        species.append(
            {
                "id": species_id,
                "chromosomes": [
                    {
                        "id": "chr1",
                        "length": n_blocks * 100,
                        "blocks": blocks,
                        "display_rank": 1,
                        "display_orientation": (
                            -1 if index % 2 else 1
                        ),
                    }
                ],
            }
        )

    return fixture_from_dict(
        {
            "fixture_version": 1,
            "id": f"orientation_chain_{n_species}_{n_blocks}",
            "title": "Orientation-chain residual test",
            "purpose": (
                "Exercise exact residual factor elimination on a "
                "chain-structured orientation problem."
            ),
            "species": species,
            "orientation_constraints": [],
            "expected": {},
        }
    )


class ResidualFactorSolverTests(unittest.TestCase):
    def test_matches_legacy_exact_dp_on_small_fixtures(self) -> None:
        names = (
            "layout_noise_only.json",
            "forbidden_subchromosomal_pretty_solution.json",
            "fusion_chain_tree.json",
            "fusion_chain_closed_cycle.json",
            "balanced_orientation_cycle.json",
        )

        for name in names:
            fixture = load_fixture(FIXTURES / name)
            legacy = exact_optimize_layer_dp(fixture)
            residual = exact_optimize_residual_factor_graph(
                fixture,
                work_cap_per_component=250_000,
            )
            self.assertEqual(
                residual.layout.optimized_score.crossings,
                legacy.layout.optimized_score.crossings,
                name,
            )
            self.assertEqual(
                residual.layout.optimality_status,
                "proven optimum",
                name,
            )

    def test_orientation_chain_is_peeled_without_global_enumeration(self) -> None:
        fixture = orientation_chain_fixture(12, n_blocks=3)
        result = exact_optimize_residual_factor_graph(
            fixture,
            work_cap_per_component=10_000,
        )

        self.assertEqual(result.layout.optimized_score.crossings, 0)
        self.assertEqual(
            result.layout.optimality_status,
            "proven optimum",
        )
        diagnostic = result.diagnostics[0]

        # The legacy orientation enumeration would contain 2^12 complete
        # assignments. Residual elimination should instead peel the chain.
        self.assertGreaterEqual(diagnostic.leaf_eliminations, 10)
        self.assertLess(
            diagnostic.table_entries_evaluated,
            2**12,
        )
        self.assertEqual(diagnostic.min_fill_eliminations, 0)

    def test_exact_factor_tables_remove_objective_neutral_variables(self) -> None:
        fixture = orientation_chain_fixture(8, n_blocks=1)
        result = exact_optimize_residual_factor_graph(
            fixture,
            work_cap_per_component=10_000,
        )
        diagnostic = result.diagnostics[0]

        self.assertEqual(result.layout.optimized_score.crossings, 0)
        self.assertEqual(diagnostic.variable_count, 8)
        self.assertEqual(diagnostic.isolated_variable_count, 8)
        self.assertGreater(
            diagnostic.factor_scope_variables_removed,
            0,
        )

    def test_auto_solver_uses_residual_factor_elimination_first(self) -> None:
        fixture = orientation_chain_fixture(8, n_blocks=3)
        result = optimize_auto(
            fixture,
            transition_cap_per_component=10_000,
            component_workers=1,
        )

        self.assertEqual(
            result.solver,
            "exact-residual-factor-elimination",
        )
        self.assertEqual(result.layout.optimized_score.crossings, 0)


if __name__ == "__main__":
    unittest.main()
