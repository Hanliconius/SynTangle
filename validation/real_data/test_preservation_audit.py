import unittest,sys,tempfile,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmark'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from syntangle.fixtures import fixture_from_dict
from syntangle.layout import initial_layout_state,LayoutState
from visual_summary import riparian,chromosome_geometry
from preservation_audit import check_panel,require_audit
class PreservationTests(unittest.TestCase):
    def test_directed_whole_flip_passes_and_endpoint_tampering_fails(self):
        f=fixture_from_dict(dict(fixture_version=1,id='audit',title='test',purpose='test',species=[dict(id=s,chromosomes=[dict(id='1',length=100,blocks=[dict(occurrence_id=s,homology_id='h',start=10,end=40,strand='-' if s=='B' else '+')])]) for s in ('A','B')]))
        state=initial_layout_state(f);ref=state.chromosome_order['B'][0]
        for sign in (1,-1):
            st=LayoutState(state.chromosome_order,{**state.chromosome_orientation,ref:sign});fig,ax=plt.subplots()
            riparian(ax,f,st,'','',respect_block_strand=True)
            self.assertEqual(check_panel(f,st,ax,chromosome_geometry(f,st))['links'],1)
            ax.patches[0].get_path().vertices[3,0]+=.01
            with self.assertRaises(ValueError):check_panel(f,st,ax,chromosome_geometry(f,st))
            plt.close(fig)
    def test_missing_and_partial_audit_never_complete(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):require_audit(d)
            Path(d,'preservation_audit.json').write_text(json.dumps(dict(status='PASS',panels=[])))
            with self.assertRaises(ValueError):require_audit(d)
if __name__=='__main__':unittest.main()
