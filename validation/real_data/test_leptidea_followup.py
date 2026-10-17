import json,random,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from syntangle.fixtures import fixture_from_dict
from syntangle.layout import initial_layout_state,score_crossings
from preservation_audit import digest
from leptidea_followup import seed,collect,CASE
class FollowupTests(unittest.TestCase):
    def build(self,root):
        perm=list(range(60));random.Random(9).shuffle(perm)
        def species(sp,order):return dict(id=sp,chromosomes=[dict(id='chr1',length=600,display_rank=1,blocks=[dict(occurrence_id=sp+str(h),homology_id=str(h),start=i*10,end=i*10+5,strand='+') for i,h in enumerate(order)])])
        raw=dict(fixture_version=1,id=CASE,title='test',purpose='test',species=[species('A',list(range(60))),species('B',perm)])
        f=fixture_from_dict(raw);state=initial_layout_state(f);upper=score_crossings(f,state).crossings
        self.assertTrue(359<=upper<=3468)
        source=root/'source';prior=root/'prior';out=source/'prepared'/CASE;out.mkdir(parents=True);fp=out/'fixture.json';fp.write_text(json.dumps(raw));sha=digest(fp)
        d=source/'figures'/CASE;d.mkdir(parents=True);(d/'comparison_audit.json').write_text(json.dumps(dict(fixture_sha256=sha,panels=[dict(title='baseline',state=state.to_dict(),crossings=upper)]*4,provenance={})))
        d=prior/'figures'/CASE;d.mkdir(parents=True);(d/'fixture.json').write_bytes(fp.read_bytes())
        trial=dict(fixture_sha256=sha,upper_bound=upper,lower_bound=359,optimized_state=state.to_dict())
        (d/'summary.json').write_text(json.dumps(dict(fixture_sha256=sha,upper_bound=upper,lower_bound=359,trials=[trial])))
        return source,prior,trial
    @patch('leptidea_followup.require_audit')
    def test_changed_fixture_is_rejected(self,_):
        with tempfile.TemporaryDirectory() as d:
            source,prior,trial=self.build(Path(d));seed(source,prior)
            p=source/'prepared'/CASE/'fixture.json';p.write_text(p.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'different evidence'):seed(source,prior)
    @patch('leptidea_followup.require_audit')
    @patch('leptidea_followup.render')
    def test_weaker_new_bounds_cannot_replace_prior_bound(self,*_):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);source,prior,trial=self.build(root);run=root/'run';run.mkdir()
            for i in range(2):
                out=run/f'task_{i}';out.mkdir();(out/'COMPLETE').write_text('PASS');(out/'result.json').write_text(json.dumps(dict(trial,lower_bound=0)))
            collect(source,prior,run)
            final=json.loads((run/'figures'/CASE/'summary.json').read_text())
            self.assertEqual(final['lower_bound'],359);self.assertEqual(final['upper_bound'],trial['upper_bound'])
if __name__=='__main__':unittest.main()
