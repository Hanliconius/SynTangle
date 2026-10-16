import json,sys,tempfile,unittest,subprocess
from pathlib import Path
from unittest.mock import patch
from syntangle.fixtures import fixture_from_dict
from syntangle.layout import initial_layout_state
import plot_comparison as p

class MissingNativeTests(unittest.TestCase):
    def test_failed_native_still_runs_and_audits_optimizer(self):
        fixture=fixture_from_dict(dict(fixture_version=1,id='test',title='t',purpose='t',species=[dict(id=sp,chromosomes=[dict(id='1',length=100,blocks=[dict(occurrence_id=sp,homology_id='h',start=10,end=20,strand='+')])]) for sp in ('a','b')]))
        state=initial_layout_state(fixture)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'fixture.json').write_text('{}');(root/'provenance.json').write_text(json.dumps(dict(retained_links=1,pair_audit=[])))
            result=dict(optimized_state=state.to_dict(),upper_bound=0,lower_bound=0,optimality_status='test')
            with patch.object(sys,'argv',['plot','--prepared',d,'--output',str(root/'out'),'--seconds','1']),patch.object(p,'load_fixture',return_value=fixture),patch.object(p,'export_native_input'),patch.object(p,'add_deposited_reference_anchors',return_value={}),patch.object(p.subprocess,'run',side_effect=subprocess.CalledProcessError(1,'Rscript')),patch.object(p,'run_hybrid',return_value=result) as run:
                p.main()
            report=json.loads((root/'out/comparison_audit.json').read_text())
            self.assertEqual(report['native_status'],'unavailable');self.assertIsNone(report['panels'][1]['state'])
            self.assertEqual(len(p.require_audit(root/'out')['panels']),2)
            self.assertTrue((root/'out/COMPLETE').exists());run.assert_called_once()
if __name__=='__main__':unittest.main()
