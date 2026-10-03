from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import (
    build_incidence_graph,
    exact_optimize_layer_dp,
    load_fixture,
    optimize_branch_and_bound,
)
from syntangle.bounds import build_relaxed_crossing_bound
from syntangle.incidence import chromosome_node_id
from syntangle.model import (
    BlockOccurrence,
    Chromosome,
    ChromosomeRef,
    Fixture,
)
from syntangle.orientation_space import orientation_basis


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


def _block(
    occurrence_id: str,
    homology_id: str,
    position: float,
) -> BlockOccurrence:
    return BlockOccurrence(
        occurrence_id=occurrence_id,
        homology_id=homology_id,
        start=position,
        end=position + 1.0,
        strand="+",
    )


def high_orientation_fixture() -> Fixture:
    a_ref = ChromosomeRef("speciesA", "A0")
    homologies = [f"H{i}" for i in range(1, 15)]

    a = Chromosome(
        ref=a_ref,
        length=200.0,
        display_rank=0,
        display_orientation=1,
        blocks=tuple(
            _block(f"A_{homology}", homology, i * 10.0)
            for i, homology in enumerate(homologies, start=1)
        ),
    )

    b_chromosomes = []
    b1_ref = ChromosomeRef("speciesB", "B1")
    b_chromosomes.append(
        Chromosome(
            ref=b1_ref,
            length=100.0,
            display_rank=0,
            display_orientation=1,
            blocks=(
                _block("B_H1", "H1", 10.0),
                _block("B_H3", "H3", 20.0),
                _block("B_H2", "H2", 30.0),
            ),
        )
    )

    for index, homology in enumerate(homologies[3:], start=2):
        ref = ChromosomeRef("speciesB", f"B{index}")
        b_chromosomes.append(
            Chromosome(
                ref=ref,
                length=100.0,
                display_rank=index - 1,
                display_orientation=1,
                blocks=(
                    _block(
                        f"{ref.chromosome_id}_{homology}",
                        homology,
                        50.0,
                    ),
                ),
            )
        )

    return Fixture(
        fixture_version=1,
        fixture_id="high_orientation_monotone",
        title="High free-orientation regression",
        purpose=(
            "Ensure compact orientation bases are not expanded before "
            "branch-and-bound"
        ),
        chromosomes=(a, *b_chromosomes),
        orientation_constraints=(),
        expected={},
    )


class MonotoneReductionTests(unittest.TestCase):
    def test_relaxed_bound_never_exceeds_exact_optimum(self) -> None:
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
            graph = build_incidence_graph(fixture)
            total_lower = 0

            for component in graph.connected_components():
                refs = tuple(
                    sorted(
                        ref
                        for ref in fixture.chromosome_refs
                        if chromosome_node_id(ref) in component
                    )
                )
                basis = orientation_basis(fixture, refs)
                bound = build_relaxed_crossing_bound(
                    fixture,
                    component,
                    basis,
                )
                total_lower += bound.lower_bound(
                    tuple(None for _ in basis.free_flip_groups)
                )

            self.assertLessEqual(
                total_lower,
                exact.layout.optimized_score.crossings,
                name,
            )

    def test_more_than_4096_orientations_remain_compact(self) -> None:
        fixture = high_orientation_fixture()
        graph = build_incidence_graph(fixture)
        components = graph.connected_components()
        self.assertEqual(len(components), 1)

        refs = tuple(sorted(fixture.chromosome_refs))
        basis = orientation_basis(fixture, refs)
        self.assertEqual(len(basis.free_flip_groups), 13)
        self.assertEqual(basis.assignment_count, 8192)

        result = optimize_branch_and_bound(
            fixture,
            node_cap_per_component=50,
            orientation_cap_per_component=1,
            local_restarts=2,
            seed=53,
        )

        self.assertEqual(result.upper_bound, 1)
        self.assertEqual(result.lower_bound, 1)
        self.assertEqual(result.gap, 0)
        self.assertEqual(
            result.layout.optimality_status,
            "proven optimum",
        )

        diagnostic = result.components[0]
        self.assertEqual(diagnostic.implicit_orientation_states, 8192)
        self.assertEqual(diagnostic.root_lower_bound, 1)
        self.assertEqual(diagnostic.orientation_assignments, 0)
        self.assertEqual(diagnostic.orientation_nodes_evaluated, 0)

    def test_node_cutoff_keeps_valid_global_bound(self) -> None:
        fixture = load_fixture(
            FIXTURES / "fusion_chain_closed_cycle.json"
        )
        exact = exact_optimize_layer_dp(fixture)
        result = optimize_branch_and_bound(
            fixture,
            node_cap_per_component=1,
            local_restarts=2,
            seed=59,
        )

        optimum = exact.layout.optimized_score.crossings
        self.assertLessEqual(result.lower_bound, optimum)
        self.assertGreaterEqual(result.upper_bound, optimum)


if __name__ == "__main__":
    unittest.main()
