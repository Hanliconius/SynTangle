from pathlib import Path
import unittest
from unittest.mock import patch

from syntangle import load_fixture, optimize_auto, score_crossings
from syntangle.heuristic import _component_map
from syntangle.layout import initial_layout_state
from syntangle.search_control import _progress, _deadline
import syntangle.branch_bound as branch

ROOT = Path(__file__).resolve().parents[1]


class DeadlineTests(unittest.TestCase):
    def setUp(self):
        self.fixture = load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')

    def test_unstarted_components_retain_incumbent_and_valid_bounds(self):
        with patch('syntangle.search_pipeline.expired', return_value=True):
            result = optimize_auto(self.fixture, local_restarts=1, time_limit_seconds=1)
        self.assertEqual(result.details['lower_bound'], 0)
        self.assertEqual(result.details['upper_bound'],
                         score_crossings(self.fixture, result.layout.optimized_state).crossings)
        self.assertTrue(all(d['method'] == 'deadline-unstarted'
                            for d in result.details['component_diagnostics']))
        self.assertIsNone(_deadline.get())
        self.assertIsNone(_progress.get())

    def test_improvements_stream_from_inside_branch_search(self):
        components, component_of = _component_map(self.fixture)
        state = initial_layout_state(self.fixture)
        seen = []
        token = _progress.set(lambda candidate: seen.append(
            score_crossings(self.fixture, candidate, restrict_component_nodes=components[0]).crossings))
        try:
            result = branch._solve_branch_component(self.fixture, state, state, 0,
                                                    components[0], component_of, 1000)
        finally:
            _progress.reset(token)
        self.assertTrue(seen)
        self.assertEqual(seen, sorted(set(seen), reverse=True))
        self.assertEqual(seen[-1], result.diagnostic.upper_bound)

    def test_mid_search_deadline_does_not_claim_false_optimum(self):
        components, component_of = _component_map(self.fixture)
        state = initial_layout_state(self.fixture)
        exact = branch._solve_branch_component(self.fixture, state, state, 0,
                                               components[0], component_of, 10000)
        calls = 0
        def cutoff():
            nonlocal calls
            calls += 1
            return calls > 30
        with patch('syntangle.branch_bound.expired', side_effect=cutoff):
            partial = branch._solve_branch_component(self.fixture, state, state, 0,
                                                     components[0], component_of, 10000)
        self.assertLessEqual(partial.lower_bound, exact.diagnostic.upper_bound)
        self.assertGreaterEqual(partial.diagnostic.upper_bound, exact.diagnostic.upper_bound)
        self.assertEqual(partial.diagnostic.proven,
                         partial.lower_bound == partial.diagnostic.upper_bound)

    def test_long_deadline_preserves_search_and_restores_context(self):
        kwargs = dict(local_restarts=1, transition_cap_per_component=1,
                      branch_node_cap_per_component=50)
        previous = optimize_auto(self.fixture, **kwargs)
        callbacks = []
        result = optimize_auto(self.fixture, **kwargs, time_limit_seconds=100,
            progress_callback=lambda state, score: callbacks.append((state, score)))
        self.assertEqual(result.layout, previous.layout)
        self.assertEqual(result.details['component_diagnostics'], previous.details['component_diagnostics'])
        self.assertTrue(callbacks)
        for state, score in callbacks:
            self.assertEqual(score_crossings(self.fixture, state).crossings, score)
        self.assertEqual([score for _, score in callbacks],
                         sorted(set(score for _, score in callbacks), reverse=True))
        self.assertIsNone(_deadline.get())
        self.assertIsNone(_progress.get())

    def test_parallel_deadline_rejected_explicitly(self):
        with self.assertRaisesRegex(ValueError, 'component_workers=1'):
            optimize_auto(self.fixture, time_limit_seconds=1, component_workers=2)
