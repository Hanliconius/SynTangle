from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import derive_ordering_constraints, load_fixture


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "examples" / "fixtures"


class OrderingConstraintTests(unittest.TestCase):
    def test_forbidden_pretty_solution_preserves_biological_chains(self) -> None:
        fixture = load_fixture(FIXTURES / "forbidden_subchromosomal_pretty_solution.json")
        state = derive_ordering_constraints(fixture)
        chains = {chain.chromosome.label: chain for chain in state.reversible_chains}

        self.assertEqual(chains["sp1:X"].homology_ids, ("HA", "HB", "HC", "HD"))
        self.assertEqual(chains["sp1:X"].reversed_homology_ids, ("HD", "HC", "HB", "HA"))
        self.assertEqual(chains["sp2:Y"].homology_ids, ("HA", "HC", "HB", "HD"))
        self.assertEqual(chains["sp2:Y"].reversed_homology_ids, ("HD", "HB", "HC", "HA"))

        forbidden = tuple(fixture.expected["forbidden_internal_homology_orders"]["sp1:X"][0])
        self.assertNotEqual(chains["sp1:X"].homology_ids, forbidden)
        self.assertNotEqual(chains["sp1:X"].reversed_homology_ids, forbidden)

    def test_layout_noise_components_are_equivalent_not_unresolved(self) -> None:
        fixture = load_fixture(FIXTURES / "layout_noise_only.json")
        state = derive_ordering_constraints(fixture)
        self.assertEqual(len(state.equivalent_components), 4)
        self.assertEqual(len(state.unresolved_chromosome_orders), 0)

    def test_fusion_tree_exposes_within_component_unresolved_orders(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_tree.json")
        state = derive_ordering_constraints(fixture)
        unresolved = {
            (item.species_id, item.a.chromosome_id, item.b.chromosome_id)
            for item in state.unresolved_chromosome_orders
        }
        self.assertIn(("sp2", "B", "C"), unresolved)
        self.assertEqual(state.forced_internal_chain_count, 3)

    def test_no_cross_species_adjacency_is_invented(self) -> None:
        fixture = load_fixture(FIXTURES / "fusion_chain_closed_cycle.json")
        state = derive_ordering_constraints(fixture)
        # The state contains only observed within-chromosome reversible chains,
        # disconnected-component equivalence, and unresolved whole-chromosome
        # pairs. There is no inferred hard adjacency field by design.
        self.assertFalse(hasattr(state, "forced_cross_species_adjacencies"))


if __name__ == "__main__":
    unittest.main()
