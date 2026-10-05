import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest

from syntangle import load_fixture, optimize_auto

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('pipeline_diagnostic',
    ROOT/'validation/benchmark/diagnose_pipeline.py')
diagnostic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(diagnostic)


class DiagnosticTests(unittest.TestCase):
    def test_instrumentation_preserves_solver_result_and_restores_functions(self):
        import syntangle.search_pipeline as pipeline
        original = pipeline.median_order_seed
        fixture = load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        kwargs = dict(transition_cap_per_component=1, branch_node_cap_per_component=1000,
                      local_restarts=1, local_max_improving_steps=5, seed=32)
        expected = optimize_auto(fixture, **kwargs).to_dict()
        observer = diagnostic.Observer().install()
        try:
            observed = optimize_auto(fixture, **kwargs).to_dict()
        finally:
            observer.close()
        self.assertEqual(expected, observed)
        self.assertIs(pipeline.median_order_seed, original)
        self.assertTrue(observer.stats['median_order_seed']['completed'])
        self.assertTrue(observer.stats['_solve_branch_component']['completed'])
        self.assertTrue(any(s['event'].endswith('.exception') and
                            s['type'] == 'SearchSpaceTooLarge' for s in observer.samples))

    def test_orientation_reduction_tuple_return_preserved_and_counted(self):
        import syntangle.branch_bound as branch
        class Bound:
            def lower_bound(self, bits):
                # Bit zero is provably unable to improve the incumbent.
                return 5 if bits[0] == 0 else 0
        kwargs = dict(incumbent_upper=5, bound=Bound(), cache={})
        expected = branch._reduce_orientation((None,), **kwargs)
        observer = diagnostic.Observer().install()
        try:
            actual = branch._reduce_orientation((None,), incumbent_upper=5,
                                                bound=Bound(), cache={})
        finally:
            observer.close()
        self.assertEqual(actual, expected)
        reduced, hits = actual
        self.assertEqual(reduced.bits, (1,))
        self.assertEqual(reduced.forced, 1)
        self.assertEqual(observer.stats['_reduce_orientation']['forced_total'], 1)
        self.assertEqual(observer.stats['_reduce_orientation']['cache_hits'], hits)
        self.assertEqual(observer.samples[-1]['event'], 'orientation.reduction')

    def test_audit_distinguishes_unresolved_timeout_from_completed_bounds(self):
        fixture = load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        result = optimize_auto(fixture, local_restarts=1).to_dict()
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)/'pipeline_pair_unit'/'case'
            case.mkdir(parents=True)
            (case/'paired.json').write_text(json.dumps(dict(case_id='complete', methods=[
                dict(method='pipeline', status='complete', result=result)])))
            timeout = case.parent/'timeout'
            timeout.mkdir()
            (timeout/'paired.json').write_text(json.dumps(dict(case_id='timeout', methods=[
                dict(method='pipeline', status='timeout', incumbent=dict(crossings=4, seconds=1))])))
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                diagnostic.audit(directory)
            self.assertIn('bounds=0..0 gap=0', output.getvalue())
            self.assertIn('NO FINAL REDUCTION AUDIT OR PROOF GAP', output.getvalue())


if __name__ == '__main__':
    unittest.main()
