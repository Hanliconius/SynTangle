import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmark'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from syntangle.fixtures import fixture_from_dict
from syntangle.layout import initial_layout_state,LayoutState
from visual_summary import riparian

class DirectedRibbonTests(unittest.TestCase):
    def test_whole_flip_reflects_each_directed_endpoint(self):
        for strand in ('+','-'):
            fixture=fixture_from_dict(dict(fixture_version=1,id='directed_flip',title='test',purpose='test',species=[dict(id=sp,chromosomes=[dict(id='1',length=100,blocks=[dict(occurrence_id=sp,homology_id='h',start=10,end=30,strand=strand)])]) for sp in ('A','B')]))
            state=initial_layout_state(fixture);ref=state.chromosome_order['B'][0]
            flipped=LayoutState(state.chromosome_order,{**state.chromosome_orientation,ref:-1})
            fig,axes=plt.subplots(1,2)
            riparian(axes[0],fixture,state,'','',respect_block_strand=True)
            riparian(axes[1],fixture,flipped,'','',respect_block_strand=True)
            v=axes[0].patches[0].get_path().vertices;w=axes[1].patches[0].get_path().vertices
            # One chromosome per row: normalized extent 0..1. Reflection
            # preserves endpoint identity, not just the unordered interval.
            self.assertAlmostEqual(w[3,0],1-v[3,0])
            self.assertAlmostEqual(w[4,0],1-v[4,0])
            self.assertEqual(fixture.chromosomes[1].blocks[0].start,10)
            self.assertEqual(fixture.chromosomes[1].blocks[0].end,30)
            plt.close(fig)
if __name__=='__main__':unittest.main()
