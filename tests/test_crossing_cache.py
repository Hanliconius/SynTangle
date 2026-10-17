import random
from dataclasses import replace
import unittest
from pathlib import Path
from unittest.mock import patch
from syntangle import load_fixture, initial_layout_state, score_crossings
from syntangle.crossing_cache import CrossingCostCache
from syntangle.layout import LayoutState, AmbiguousHomologyError
from syntangle.heuristic import optimize_local_search

ROOT = Path(__file__).resolve().parents[1]

class ReferenceCache:
    def __init__(self, fixture):
        self.fixture = fixture
    def prepare(self, state):
        return score_crossings(self.fixture, state).crossings
    score_candidate = prepare

class CachedCrossingTests(unittest.TestCase):
    def test_random_orders_and_multi_chromosome_flips_match_full_scorer(self):
        for path in sorted((ROOT / 'examples/fixtures').glob('*.json')):
            if path.name == 'fixture.schema.json':
                continue
            fixture = load_fixture(path)
            try:
                cache = CrossingCostCache(fixture)
            except AmbiguousHomologyError:
                continue  # Some validation fixtures intentionally contain ambiguous homology.
            rng = random.Random(48)
            state = initial_layout_state(fixture)
            self.assertEqual(cache.prepare(state), score_crossings(fixture, state).crossings)
            for _ in range(40):
                orders = {}
                for species, refs in state.chromosome_order.items():
                    refs = list(refs)
                    rng.shuffle(refs)
                    orders[species] = tuple(refs)
                candidate = LayoutState(orders, {ref: rng.choice((-1, 1)) for ref in state.chromosome_orientation})
                expected = score_crossings(fixture, candidate).crossings
                self.assertEqual(cache.score_candidate(candidate), expected, path.name)
                self.assertEqual(cache.prepare(candidate), expected, path.name)
                state = candidate

    def test_identical_search_decisions_and_evaluation_counts(self):
        for name in ('fusion_chain_tree.json', 'fusion_chain_closed_cycle.json'):
            path = ROOT / 'examples/fixtures' / name
            if not path.exists():
                continue
            if path.name == 'fixture.schema.json':
                continue
            fixture = load_fixture(path)
            for seed in range(3):
                fast = optimize_local_search(fixture, restarts=2, max_improving_steps=5, seed=seed)
                with patch('syntangle.heuristic.CrossingCostCache', ReferenceCache):
                    reference = optimize_local_search(fixture, restarts=2, max_improving_steps=5, seed=seed)
                self.assertEqual(fast.to_dict(), reference.to_dict())

    def test_tied_positions_and_reflection_rounding(self):
        fixture = load_fixture(ROOT / 'examples/fixtures/perfect_1to1_30x3.json')
        chromosomes = []
        for chromosome in fixture.chromosomes:
            blocks = tuple(replace(block, start=1e-20 if i % 2 else 0.0,
                                   end=1e-20 if i % 2 else 0.0)
                           for i, block in enumerate(chromosome.blocks))
            chromosomes.append(replace(chromosome, blocks=blocks))
        fixture = replace(fixture, chromosomes=tuple(chromosomes))
        cache = CrossingCostCache(fixture)
        state = initial_layout_state(fixture)
        cache.prepare(state)
        for sign in (-1, 1):
            candidate = LayoutState(state.chromosome_order,
                {ref: sign for ref in state.chromosome_orientation})
            self.assertEqual(cache.score_candidate(candidate), score_crossings(fixture, candidate).crossings)

    def test_progress_scores_are_valid_and_strictly_improve(self):
        fixture = load_fixture(ROOT / 'examples/fixtures/fusion_chain_tree.json')
        scores = []
        def report(state, score):
            self.assertEqual(score, score_crossings(fixture, state).crossings)
            scores.append(score)
        optimize_local_search(fixture, restarts=2, seed=4, progress_callback=report)
        self.assertTrue(scores)
        self.assertTrue(all(a > b for a, b in zip(scores, scores[1:])))
