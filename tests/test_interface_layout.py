import unittest
from pathlib import Path
from syntangle.table_input import load_tables
from syntangle.layout import initial_layout_state, score_crossings
from syntangle.interface_layout import optimize_interface
from syntangle.saved_layout import decode_saved_layout

class Tests(unittest.TestCase):
    def fixture(self):
        p=Path(__file__).resolve().parents[1]/'examples/quickstart'
        return load_tables(p/'chromosomes.tsv',p/'links.tsv')[0]
    def test_transaction_and_full_evidence(self):
        f=self.fixture(); start=initial_layout_state(f); calls=[]
        def solver(sub,seed,*args,**kwargs):
            calls.append(sub.species_ids)
            return {'optimized_state':seed.to_dict(),'status':'mock unchanged, no certificate'}
        result=optimize_interface(f,start,seconds=2,row_seconds=.05,pair_solver=solver)
        end=decode_saved_layout(f,result['optimized_state'])
        self.assertEqual(len(calls[0]),2)
        self.assertLessEqual(score_crossings(f,end).crossings,score_crossings(f,start).crossings)
        self.assertEqual(set(end.chromosome_orientation),set(f.chromosome_refs))
        self.assertEqual(result['input']['link_count'],result['output']['link_count'])
        self.assertIsNone(result['global_lower_bound'])
        self.assertEqual(result['accepted'],result['proposal_score']['crossings']<result['input']['crossings'])
    def test_full_score_rejects_locally_good_pair(self):
        from syntangle.layout import LayoutState, canonicalize_component_order
        import time
        f=self.fixture(); start=initial_layout_state(f)
        order={sp:tuple(sorted(refs)) for sp,refs in start.chromosome_order.items()}
        start=canonicalize_component_order(f,LayoutState(order,start.chromosome_orientation))
        def solver(sub,seed,*args,**kwargs):
            changed=LayoutState({sp:tuple(reversed(refs)) for sp,refs in seed.chromosome_order.items()},
                                {ref:-sign for ref,sign in seed.chromosome_orientation.items()})
            time.sleep(.003)
            return {'optimized_state':changed.to_dict()}
        result=optimize_interface(f,start,seconds=.001,boundary=('A','B'),pair_solver=solver)
        self.assertFalse(result['accepted'])
        self.assertGreater(result['proposal_score']['crossings'],result['output']['crossings'])
    def test_invalid_budget(self):
        f=self.fixture()
        with self.assertRaises(ValueError):optimize_interface(f,initial_layout_state(f),seconds=-1)
if __name__=='__main__':unittest.main()
