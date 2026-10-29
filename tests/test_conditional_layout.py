import itertools,json,tempfile,unittest
from pathlib import Path
from dataclasses import replace
from syntangle.table_input import load_tables
from syntangle.layout import initial_layout_state,score_crossings,LayoutState,canonicalize_component_order
from syntangle.heuristic import _component_map
from syntangle.conditional_layout import PairScorer,align_blocks,mirror_block,row_candidates,optimize_conditional
from syntangle.orientation_space import orientation_basis

E=Path(__file__).resolve().parents[1]/'examples/quickstart'
class Tests(unittest.TestCase):
 def fixture(self):return load_tables(E/'chromosomes.tsv',E/'links.tsv')[0]
 def test_pair_delta_matches_full_objective(self):
  f=self.fixture();state=canonicalize_component_order(f,initial_layout_state(f));s=PairScorer(f);_,comp=_component_map(f);basis=orientation_basis(f,f.chromosome_refs)
  for sp in f.species_ids:
   for c in row_candidates(state,sp,comp,basis):
    self.assertEqual(s.score(c,[sp])-s.score(state,[sp]),score_crossings(f,c).crossings-score_crossings(f,state).crossings)
    for other in f.species_ids:
     if other!=sp:self.assertEqual(c.chromosome_order[other],state.chromosome_order[other])
  self.assertEqual(s.score(state),score_crossings(f,state).crossings)
 def test_binary_block_dp_matches_exhaustive(self):
  f=self.fixture();state=canonicalize_component_order(f,initial_layout_state(f));blocks=[[sp] for sp in f.species_ids];_,comp=_component_map(f);scorer=PairScorer(f)
  final,info=align_blocks(f,state,blocks,scorer,comp);scores=[]
  for bits in itertools.product((0,1),repeat=len(blocks)):
   st=state
   for b,bit in zip(blocks,bits):
    if bit:st=mirror_block(st,b,comp)
   scores.append(score_crossings(f,st).crossings)
  self.assertEqual(score_crossings(f,final).crossings,min(scores))
  with self.assertRaises(ValueError):align_blocks(f,state,list(reversed(blocks)),scorer,comp)
 def test_complete_pipeline_has_live_events_and_preserves_evidence(self):
  f=self.fixture();start=initial_layout_state(f);events=[];points=[];original=f
  result=optimize_conditional(f,start,[[sp] for sp in f.species_ids],seconds=5,row_seconds=.3,sweeps=3,event=events.append,progress=points.append)
  self.assertLess(result['output']['crossings'],result['input']['crossings']);self.assertEqual(result['input']['link_count'],result['output']['link_count'])
  self.assertIs(f,original);self.assertIsNone(result['global_lower_bound']);self.assertIn('heuristic',result['status'])
  self.assertTrue(events);self.assertTrue(points)
  for p in points:self.assertEqual(set(p.chromosome_orientation),set(f.chromosome_refs))
 def test_no_global_model_is_constructed(self):
  from unittest.mock import patch
  f=self.fixture()
  with patch('syntangle.global_experiment.JointModel',side_effect=AssertionError('global model built')):
   optimize_conditional(f,initial_layout_state(f),[[sp] for sp in f.species_ids],seconds=.5,row_seconds=.05)
 def test_bad_budget_and_blocks(self):
  f=self.fixture()
  with self.assertRaises(ValueError):optimize_conditional(f,initial_layout_state(f),[[sp] for sp in f.species_ids],seconds=float('nan'))
if __name__=='__main__':unittest.main()
