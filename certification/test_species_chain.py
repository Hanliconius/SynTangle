import itertools
import random
import unittest
from proof import Evidence
from species_chain import construct, construct_with_separators
from structural import lower_bound


def random_chain(sizes, seed=1):
    rng = random.Random(seed)
    rows = [dict(id=f'S{i}',chromosomes=[dict(id=str(j),length=100,blocks=[])
                                      for j in range(size)]) for i, size in enumerate(sizes)]
    for i in range(len(rows)-1):
        for h in range(4):
            homology = f'{i}_{h}'
            for row in (rows[i], rows[i+1]):
                c = rng.choice(row['chromosomes'])
                x = rng.randrange(40)
                c['blocks'].append(dict(occurrence_id=row['id']+homology,homology_id=homology,
                                        start=2*x,end=2*x+1,strand='+'))
    return dict(fixture_version=1,id='chain',species=rows)


def conflicting_middle(four=False):
    rows = [dict(id='A',chromosomes=[dict(id='hub',length=100,blocks=[])]),
            dict(id='B',chromosomes=[dict(id=c,length=100,blocks=[]) for c in 'abc']),
            dict(id='C',chromosomes=[dict(id='hub',length=100,blocks=[])])]
    for side, partner_order in [(0,'abc'),(2,'acb')]:
        for rank, label in enumerate(partner_order):
            h = str(side)+label
            rows[side]['chromosomes'][0]['blocks'].append(dict(occurrence_id='hub'+h,homology_id=h,
                start=2*rank,end=2*rank+1,strand='+'))
            rows[1]['chromosomes']['abc'.index(label)]['blocks'].append(dict(occurrence_id='leaf'+h,
                homology_id=h,start=10+side*10,end=11+side*10,strand='+'))
    if four:
        block = dict(occurrence_id='C_cd',homology_id='cd',start=50,end=51,strand='+')
        rows[2]['chromosomes'][0]['blocks'].append(block)
        rows.append(dict(id='D',chromosomes=[dict(id='hub',length=100,
            blocks=[{**block,'occurrence_id':'D_cd'}])]))
    return dict(fixture_version=1,id='shared_middle_conflict',species=rows)


def orientation_conflict():
    raw = dict(fixture_version=1,id='two_partner_middle_conflict',species=[
        dict(id='A',chromosomes=[dict(id='hub',length=100,blocks=[])]),
        dict(id='B',chromosomes=[dict(id=c,length=100,blocks=[]) for c in 'ab']),
        dict(id='C',chromosomes=[dict(id='hub',length=100,blocks=[])])])
    for side in (0,2):
        for partner,label in enumerate('ab'):
            for j in range(2):
                h=f'{side}_{label}_{j}'
                slot=1-j if side==2 and label=='b' else j
                x=10*partner+2*slot
                raw['species'][side]['chromosomes'][0]['blocks'].append(dict(
                    occurrence_id='hub'+h,homology_id=h,start=x,end=x+1,strand='+'))
                x=10*side+2*j
                raw['species'][1]['chromosomes'][partner]['blocks'].append(dict(
                    occurrence_id='leaf'+h,homology_id=h,start=x,end=x+1,strand='+'))
    return raw


def wide_endpoint_chain():
    from test_fusion_theorem import star
    positions={'a':[2,4,9],'b':[1,6,8],'c':[3,5,7]}
    positions.update({f'z{i:03}':[10+i] for i in range(97)})
    rows={sp:dict(id=sp,chromosomes=[]) for sp in 'ABCD'}
    for prefix,mapping in [('left',{'A':'B','B':'A'}),('right',{'A':'C','B':'D'})]:
        for row in star(positions)['species']:
            sp=mapping[row['id']]
            for chromosome in row['chromosomes']:
                chromosome['length']=1000
                for block in chromosome['blocks']:
                    block['occurrence_id']=prefix+block['occurrence_id']
                    block['homology_id']=prefix+block['homology_id']
                rows[sp]['chromosomes'].append(chromosome)
    for i in range(3):
        for sp,slot in [('B',i),('C',[1,2,0][i])]:
            rows[sp]['chromosomes'][0]['blocks'].append(dict(occurrence_id=sp+str(i),homology_id='middle'+str(i),
                start=400+10*slot,end=401+10*slot,strand='+'))
    return dict(fixture_version=1,id='wide_endpoints',species=list(rows.values()))


class SpeciesChainTests(unittest.TestCase):
    def exhaustive(self, raw):
        evidence = Evidence(raw)
        result = construct(evidence)
        self.assertEqual(result['status'], 'EXACT_SPECIES_CHAIN_RECURRENCE')
        self.assertEqual(result['upper'], min(evidence.score(s) for s in evidence.states()))
        self.assertEqual(result['lower'], result['upper'])
        return result

    def test_three_and_four_species_arbitrary_connections(self):
        for sizes, seed in itertools.product(([2,2,2],[2,2,2,2]), range(4)):
            with self.subTest(sizes=sizes, seed=seed):
                self.exhaustive(random_chain(sizes, seed))

    def test_pairwise_optima_do_not_compose(self):
        for four in (False, True):
            result = self.exhaustive(conflicting_middle(four))
            self.assertEqual(result['interface_minimum_sum'], 0)
            self.assertEqual(result['upper'], 1)
            self.assertEqual(result['compatibility_penalty'], 1)

    def test_two_middle_chromosomes_already_allow_orientation_conflict(self):
        result = self.exhaustive(orientation_conflict())
        self.assertEqual(result['interface_minimum_sum'],0)
        self.assertEqual(result['upper'],1)
        self.assertEqual(lower_bound(Evidence(orientation_conflict()))['lower'],1)

    def test_singleton_separator_gluing_matches_full_enumeration(self):
        for sizes in ([2,1,2], [2,1,2,2], [2,1,1,2]):
            for seed in range(3):
                raw=random_chain(sizes,seed)
                result=self.exhaustive(raw)
                split=construct_with_separators(raw)
                self.assertEqual(split['upper'],result['upper'])
                self.assertEqual(split['lower'],result['lower'])

    def test_long_singleton_chain_and_nonlocal_equations(self):
        raw=random_chain([1]*50)
        result=construct_with_separators(raw)
        self.assertEqual(result['upper'],construct(Evidence(raw))['upper'])
        self.assertEqual(len(result['pieces']),49)
        raw['orientation_constraints']=[dict(a='S0:0',b='S49:0',xor=0)]
        self.assertEqual(construct_with_separators(raw)['status'],'OUTSIDE_LOCAL_ORIENTATION_CONSTRAINTS')

    def test_four_species_202_chromosomes_glue_three_exact_pieces(self):
        raw=wide_endpoint_chain()
        self.assertEqual(construct(Evidence(raw))['status'],'UNRESOLVED_ROW_STATE_LIMIT')
        result=construct_with_separators(raw,max_cyclic_size=3)
        self.assertEqual(result['upper'],27)
        self.assertEqual(result['lower'],27)
        self.assertEqual([p['result']['upper'] for p in result['pieces']],[13,1,13])
        self.assertEqual(len(result['layout']['chromosome_orientation']),202)
        self.assertTrue(any(p['mirrored'] for p in result['pieces']))

    def test_fifty_species_uses_linear_number_of_transitions(self):
        raw = random_chain([1]*50)
        evidence = Evidence(raw)
        result = construct(evidence)
        self.assertEqual(result['upper'], lower_bound(evidence)['base'])
        self.assertEqual(result['row_state_counts'], [2]*50)
        self.assertEqual(result['transitions'], 49*4)

    def test_local_hard_equations_and_infeasibility(self):
        raw = random_chain([2,2,2])
        raw['orientation_constraints'] = [dict(a='S0:0',b='S0:1',xor=1),
                                         dict(a='S0:1',b='S1:0',xor=0)]
        self.exhaustive(raw)
        raw['orientation_constraints'].append(dict(a='S0:0',b='S0:1',xor=0))
        self.assertEqual(construct(Evidence(raw))['status'], 'INFEASIBLE_LOCAL_EQUATIONS')

    def test_limits_and_nonlocal_equations_are_explicit(self):
        raw = random_chain([2,2,2])
        self.assertEqual(construct(Evidence(raw),7)['status'], 'UNRESOLVED_ROW_STATE_LIMIT')
        self.assertEqual(construct(Evidence(raw),max_transitions=1)['status'], 'UNRESOLVED_TRANSITION_LIMIT')
        raw['orientation_constraints'] = [dict(a='S0:0',b='S2:0',xor=0)]
        self.assertEqual(construct(Evidence(raw))['status'], 'OUTSIDE_LOCAL_ORIENTATION_CONSTRAINTS')


if __name__ == '__main__':
    unittest.main()
