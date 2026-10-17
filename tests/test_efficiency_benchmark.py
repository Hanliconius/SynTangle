import importlib.util
from pathlib import Path
import tempfile
import unittest

path=Path(__file__).resolve().parents[1]/'validation/benchmark/benchmark_efficiency.py'
import sys
sys.path.insert(0,str(path.parent))
spec=importlib.util.spec_from_file_location('efficiency_benchmark',path)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class EfficiencyDesignTests(unittest.TestCase):
    def test_factorial_design_balanced_and_no_duplicate_tasks(self):
        entries=[dict(case_id=f'case_{i}') for i in range(9)]
        tasks=m.design(entries)
        self.assertEqual(len(tasks),252)
        keys={(r['case_id'],r['variant'],r['start'],r['budget'],r['repeat']) for r in tasks}
        self.assertEqual(len(keys),252)
        self.assertEqual(sum(r['budget']==150 for r in tasks),216)
        self.assertEqual(sum(r['budget']==600 for r in tasks),36)
        self.assertTrue(all(r['case_index']>=6 for r in tasks if r['budget']==600))

    def test_input_fingerprint_detects_changed_evidence(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);(root/'input.tsv').write_text('a\n')
            before=m.fingerprint(root);(root/'input.tsv').write_text('b\n')
            self.assertNotEqual(before,m.fingerprint(root))

class HybridAuditTests(unittest.TestCase):
    def test_new_hybrid_optimum_is_found_and_audited(self):
        import json
        from unittest.mock import patch
        from test_coupled_bound import fixture_for
        from syntangle.global_experiment import JointModel
        from syntangle.layout import initial_layout_state,LayoutState,score_crossings
        from syntangle.incidence import chromosome_node_id
        from audit_global_results import audit
        f=fixture_for(19)
        model=JointModel(f,frozenset(chromosome_node_id(r) for r in f.chromosome_refs),f.chromosome_refs)
        raw=initial_layout_state(f);start=LayoutState(raw.chromosome_order,model.basis.base_assignment)
        optimum,bound,_=model.highs(start,5)
        score=score_crossings(f,optimum).crossings
        self.assertEqual(score,bound)
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);(root/'benchmark_manifest.tsv').write_text('case_id\tcase_dir\nnew\tnew\n')
            saved=root/'hybrid_methods_only'/'new'/'hybrid_adaptive_600';saved.mkdir(parents=True)
            (saved/'result.json').write_text(json.dumps(dict(optimized_state=optimum.to_dict(),upper_bound=score,lower_bound=bound)))
            with patch('audit_global_results.load_validation_bundle',return_value=f):
                rows=audit(root,root,root/'audit.json',seconds=5,case_index=0)
            self.assertEqual(rows[0]['status'],'strict-improvement infeasibility verified')
            self.assertIn('hybrid_methods_only',rows[0]['source'])
