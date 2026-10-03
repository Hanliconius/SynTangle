from __future__ import annotations

import unittest
from pathlib import Path

from syntangle import load_fixture

from experiments.ancestral_pruning.ancestral_pruning_poc import (
    analyze_fixture,
    canonical_adjacency,
    derive_adjacency_evidence,
    parse_newick,
    topology_adjacency_lower_bound,
)


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "experiments" / "ancestral_pruning" / "fixtures"


class AncestralPruningExperimentTests(unittest.TestCase):
    def test_terminal_fusion_yields_useful_root_conditioned_bound(self) -> None:
        fixture = load_fixture(FIXTURES / "terminal_fusion.json")
        tree = parse_newick("((sp1,sp2),(sp3,sp4));")

        evidence = derive_adjacency_evidence(fixture)
        self.assertEqual(len(evidence), 1)
        self.assertEqual(
            evidence[0].observations,
            {"sp1": 0, "sp2": 0, "sp3": 0, "sp4": 1},
        )

        analyses = analyze_fixture(fixture, tree)
        self.assertEqual(len(analyses), 1)
        result = analyses[0]
        self.assertEqual(result.adjacency, ("H1", "H2"))
        self.assertEqual(result.root_cost_absent, 1)
        self.assertEqual(result.root_cost_present, 2)
        self.assertEqual(result.optimal_root_states, (0,))

        absent = topology_adjacency_lower_bound(
            analyses,
            root_assignment={canonical_adjacency("H1", "H2"): 0},
            incumbent_adjacency_cost=1,
        )
        present = topology_adjacency_lower_bound(
            analyses,
            root_assignment={canonical_adjacency("H1", "H2"): 1},
            incumbent_adjacency_cost=1,
        )

        self.assertFalse(absent.prunable)
        self.assertEqual(absent.conditioned_lower_bound, 1)
        self.assertTrue(present.prunable)
        self.assertEqual(present.conditioned_lower_bound, 2)

    def test_ambiguous_root_is_not_pruned(self) -> None:
        fixture = load_fixture(FIXTURES / "ambiguous_root.json")
        tree = parse_newick("((sp1,sp2),(sp3,sp4));")
        analyses = analyze_fixture(fixture, tree)

        self.assertEqual(len(analyses), 1)
        result = analyses[0]
        self.assertEqual(result.root_cost_absent, 1)
        self.assertEqual(result.root_cost_present, 1)
        self.assertEqual(result.optimal_root_states, (0, 1))

        for state in (0, 1):
            bound = topology_adjacency_lower_bound(
                analyses,
                root_assignment={("H1", "H2"): state},
                incumbent_adjacency_cost=1,
            )
            self.assertFalse(bound.prunable)
            self.assertEqual(bound.conditioned_lower_bound, 1)

    def test_species_tree_changes_event_lower_bound(self) -> None:
        fixture = load_fixture(FIXTURES / "ambiguous_root.json")

        concordant = analyze_fixture(
            fixture,
            parse_newick("((sp1,sp2),(sp3,sp4));"),
        )[0]
        discordant = analyze_fixture(
            fixture,
            parse_newick("((sp1,sp3),(sp2,sp4));"),
        )[0]

        self.assertEqual(concordant.minimum_changes, 1)
        self.assertEqual(discordant.minimum_changes, 2)

    def test_newick_branch_lengths_are_ignored(self) -> None:
        tree = parse_newick(
            "((sp1:0.1,sp2:0.2):0.3,(sp3:0.4,sp4:0.5):0.6);"
        )
        self.assertEqual(
            set(tree.leaf_names()),
            {"sp1", "sp2", "sp3", "sp4"},
        )


if __name__ == "__main__":
    unittest.main()
