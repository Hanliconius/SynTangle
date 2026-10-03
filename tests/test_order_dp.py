from __future__ import annotations

import unittest
from itertools import permutations
from pathlib import Path

from syntangle import (
    build_incidence_graph,
    canonicalize_component_order,
    initial_layout_state,
    load_fixture,
    optimize_species_component_order,
    score_crossings,
)
from syntangle.layout import LayoutState


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class OrderSubsetDPTests(unittest.TestCase):
    def test_one_layer_dp_matches_brute_force(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_tree.json")
        state = canonicalize_component_order(
            fixture, initial_layout_state(fixture)
        )
        graph = build_incidence_graph(fixture)
        component = graph.connected_components()[0]

        optimized, result = optimize_species_component_order(
            fixture, state, "sp2", component
        )

        refs = state.chromosome_order["sp2"]
        brute_best = None
        for order in permutations(refs):
            orders = dict(state.chromosome_order)
            orders["sp2"] = tuple(order)
            candidate = LayoutState(
                chromosome_order=orders,
                chromosome_orientation=dict(state.chromosome_orientation),
            )
            score = score_crossings(fixture, candidate).crossings
            key = (score, tuple(ref.label for ref in order))
            if brute_best is None or key < brute_best:
                brute_best = key

        self.assertIsNotNone(brute_best)
        self.assertEqual(
            score_crossings(fixture, optimized).crossings,
            brute_best[0],
        )
        self.assertEqual(
            tuple(ref.label for ref in result.chromosome_order),
            brute_best[1],
        )

    def test_one_layer_dp_changes_whole_chromosome_order_only(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        state = canonicalize_component_order(
            fixture, initial_layout_state(fixture)
        )
        graph = build_incidence_graph(fixture)
        component = graph.connected_components()[0]

        before = {
            chrom.ref: tuple(
                (block.occurrence_id, block.start, block.end)
                for block in chrom.blocks
            )
            for chrom in fixture.chromosomes
        }
        optimized, _ = optimize_species_component_order(
            fixture, state, "sp3", component
        )
        after = {
            chrom.ref: tuple(
                (block.occurrence_id, block.start, block.end)
                for block in chrom.blocks
            )
            for chrom in fixture.chromosomes
        }

        self.assertEqual(before, after)
        self.assertEqual(
            set(optimized.chromosome_order["sp3"]),
            set(state.chromosome_order["sp3"]),
        )

    def test_subset_dp_state_count_is_not_factorial_enumeration(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        state = canonicalize_component_order(
            fixture, initial_layout_state(fixture)
        )
        graph = build_incidence_graph(fixture)
        component = graph.connected_components()[0]

        _, result = optimize_species_component_order(
            fixture, state, "sp3", component
        )
        self.assertEqual(len(result.chromosome_order), 4)
        self.assertLessEqual(result.subset_states_evaluated, 1 + 4 * (2 ** 4))


if __name__ == "__main__":
    unittest.main()
