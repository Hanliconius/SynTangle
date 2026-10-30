import unittest
from syntangle.recovery import connection_audit, recover_layout
from syntangle.crossing_cache import CrossingCostCache
from syntangle.layout import initial_layout_state, LayoutState, score_crossings
from syntangle.orientation_space import orientation_basis
from syntangle.saved_layout import validate_saved_layout
from syntangle.audit import fixture_fingerprint
from test_final_exhaustive import tiny, enumerate_legal


class RecoveryTests(unittest.TestCase):
    def test_fast_cache_and_relaxed_bound_against_exhaustive(self):
        for seed in range(3):
            f = tiny(seed)
            cache = CrossingCostCache(f,cell_aggregation=True)
            states = list(enumerate_legal(f))
            cache.prepare(states[0])
            for state in states:
                self.assertEqual(cache.score_candidate(state),score_crossings(f,state).crossings)
            lower = cache.constant+sum(min(c) for a,b,c in cache.terms)
            self.assertLessEqual(lower,min(score_crossings(f,s).crossings for s in states))

    def test_monotone_legal_recovery_and_replay(self):
        f = tiny(3)
        basis = orientation_basis(f,f.chromosome_refs)
        start = LayoutState(initial_layout_state(f).chromosome_order,basis.base_assignment)
        groups = [[(r.species_id,r.chromosome_id) for r in f.chromosome_refs if r.chromosome_id=='0']]
        fingerprint = fixture_fingerprint(f)
        checkpoints=[]
        final,result = recover_layout(f,start,groups,2,progress=checkpoints.append)
        previous = score_crossings(f,start).crossings
        refs = {r.label:r for r in f.chromosome_refs}
        replay = start
        for step, state in zip(result['accepted_moves'],checkpoints):
            validate_saved_layout(f,state,basis)
            self.assertLess(score_crossings(f,state).crossings,previous)
            previous = score_crossings(f,state).crossings
            replay = LayoutState({**replay.chromosome_order,**{sp:tuple(refs[x] for x in order) for sp,order in step['changed_orders'].items()}},
                                 {**replay.chromosome_orientation,**{refs[k]:v for k,v in step['changed_orientations'].items()}})
            self.assertEqual(replay,state)
        self.assertEqual(replay,final)
        self.assertEqual(fixture_fingerprint(f),fingerprint)
        frozen,result = recover_layout(f,start,groups,2,active_refs=set())
        self.assertEqual(frozen,start)

    def test_bridge_and_tied_partner_audit(self):
        rows=[]
        for a,b,n in [('Z','Z',5),('Z','a',1),('a','a',5)]:
            rows.extend(dict(species1='A',chromosome1=a,species2='B',chromosome2=b) for _ in range(n))
        edges,groups=connection_audit(rows,{('A','Z'):'Z',('B','Z'):'Z'})
        self.assertEqual(len(edges),3)
        self.assertTrue(all(e['graph_bridge'] for e in edges))
        self.assertEqual(len(groups),2)
        self.assertEqual(sum(e['links'] for e in edges),11)
        tied=[dict(species1='A',chromosome1='Z',species2='B',chromosome2=b) for b in ('a','b')]
        self.assertEqual(connection_audit(tied,{})[1],[])

if __name__=='__main__':unittest.main()
