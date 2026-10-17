import unittest
from syntangle.fixtures import fixture_from_dict
from syntangle.layout import initial_layout_state
from syntangle.ribbon_objective import ribbon_metrics,optimize_ribbons,twist_terms
class RibbonTests(unittest.TestCase):
 def fixture(self):
  return fixture_from_dict(dict(fixture_version=1,id='twists',title='t',purpose='t',species=[dict(id=s,chromosomes=[dict(id='1',length=100,blocks=[dict(occurrence_id=s,homology_id='h',start=10,end=80,strand='+')])]) for s in ['a','b']]))
 def test_twist_is_penalized_even_with_zero_midpoint_crossings(self):
  from syntangle.layout import LayoutState
  f=self.fixture();s=initial_layout_state(f);signs=dict(s.chromosome_orientation);signs[f.chromosome_refs[0]]=-1;s=LayoutState(s.chromosome_order,signs)
  self.assertEqual(ribbon_metrics(f,s)['twists'],1)
  for mode in ('lex','weighted'):
   result=optimize_ribbons(f,s,mode=mode,seconds=2)
   self.assertEqual(result['final_metrics']['crossings'],0);self.assertEqual(result['final_metrics']['twists'],0)
 def test_polynomial_matches_direct_metric(self):
  from syntangle.global_experiment import JointModel
  from syntangle.heuristic import _component_map
  from syntangle.layout import LayoutState
  from itertools import product
  f=self.fixture();s=initial_layout_state(f);comps,_=_component_map(f);model=JointModel(f,comps[0],f.chromosome_refs)
  for signs in product((-1,1),repeat=2):
   state=LayoutState(s.chromosome_order,dict(zip(f.chromosome_refs,signs)));bits=model.encode(state)
   for weighted in (False,True):
    c,l,q=twist_terms(f,model,weighted);value=c+l@bits+sum(v*bits[a]*bits[b] for (a,b),v in q.items())
    self.assertAlmostEqual(value,ribbon_metrics(f,state)['width_weighted_twists' if weighted else 'twists'])
 def test_modified_objectives_match_exhaustive_small_layouts(self):
  from itertools import product,permutations
  from syntangle.layout import LayoutState
  data=dict(fixture_version=1,id='trade',title='t',purpose='t',species=[])
  for sp in ('a','b'):
   chroms=[]
   for ci in (0,1):
    blocks=[dict(occurrence_id=sp+str(ci)+str(j),homology_id='h'+str(ci*2+j),start=(10+30*j if sp=='a' else 60-30*j),end=(20+30*j+5*ci if sp=='a' else 70-30*j+5*ci),strand='-' if sp=='b' and j==0 else '+') for j in (0,1)]
    chroms.append(dict(id=str(ci),length=100,blocks=blocks))
   data['species'].append(dict(id=sp,chromosomes=chroms))
  f=fixture_from_dict(data);start=initial_layout_state(f);states=[]
  for orders in product(*(list(permutations(start.chromosome_order[sp])) for sp in f.species_ids)):
   for signs in product((-1,1),repeat=4):states.append(LayoutState(dict(zip(f.species_ids,orders)),dict(zip(f.chromosome_refs,signs))))
  for mode,width in product(('lex','weighted'),(False,True)):
   key='width_weighted_twists' if width else 'twists';cap=ribbon_metrics(f,start)['crossings']
   metrics=[ribbon_metrics(f,st) for st in states]
   optimum=min(m[key] if mode=='lex' else m['crossings']+m[key] for m in metrics if mode!='lex' or m['crossings']<=cap)
   result=optimize_ribbons(f,start,mode=mode,width_weighted=width,seconds=3);m=result['final_metrics']
   self.assertAlmostEqual(m[key] if mode=='lex' else m['crossings']+m[key],optimum)
if __name__=='__main__':unittest.main()
