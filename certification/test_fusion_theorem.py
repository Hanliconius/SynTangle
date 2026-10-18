import itertools
import unittest

from proof import Evidence
from fusion_theorem import construct, _subset_dp


def star(positions, reversed_center_row=False, leaf_orders=None):
    center, leaves = [], []
    for label, coords in positions.items():
        native = list(range(len(coords))) if leaf_orders is None else leaf_orders[label]
        leaf = []
        for i, x in enumerate(coords):
            h = label + str(i)
            center.append(dict(occurrence_id='c'+h, homology_id=h, start=2*x, end=2*x+1, strand='+'))
            y = native.index(i) * 2
            leaf.append(dict(occurrence_id='l'+h, homology_id=h, start=y, end=y+1, strand='+'))
        leaves.append(dict(id=label, length=100, blocks=leaf))
    rows = [dict(id='A', chromosomes=[dict(id='center', length=100, blocks=center)]),
            dict(id='B', chromosomes=leaves)]
    if reversed_center_row:
        rows.reverse()
    return dict(fixture_version=1, id='fusion', species=rows)


class FusionTheoremTests(unittest.TestCase):
    def check_exhaustive(self, raw):
        evidence = Evidence(raw)
        result = construct(evidence)
        self.assertEqual(result['status'], 'EXACT_FUSION_STAR_CONSTRUCTION')
        optimum = min(evidence.score(s) for s in evidence.states())
        self.assertEqual(result['upper'], optimum)
        self.assertEqual(result['lower'], optimum)
        return result

    def test_interleaving_still_has_a_closed_form_positive_optimum(self):
        for reverse in (False, True):
            result = self.check_exhaustive(star({'a': [1, 3], 'b': [2, 4]}, reverse))
            self.assertEqual(result['bundle_lower'], 0)
            self.assertEqual(result['upper'], 1)
            self.assertEqual(result['components'][0]['method'], 'ACYCLIC_PAIR_PREFERENCES')

    def test_three_partner_cycle_is_a_sharp_attainment_obstruction(self):
        result = self.check_exhaustive(star({'a': [2, 4, 9], 'b': [1, 6, 8], 'c': [3, 5, 7]}))
        item = result['components'][0]
        self.assertEqual(item['pair_lower'], 12)
        self.assertEqual(result['upper'], 13)
        self.assertEqual(item['method'], 'EXACT_SUBSET_RECURRENCE')
        declined = construct(Evidence(star({'a': [2, 4, 9], 'b': [1, 6, 8], 'c': [3, 5, 7]})), 2)
        self.assertEqual(declined['status'], 'UNRESOLVED_CYCLIC_DEGREE_LIMIT')

    def test_tied_midpoints_are_not_crossings(self):
        result = self.check_exhaustive(star({'a': [1, 3], 'b': [1, 2]}))
        self.assertEqual(result['upper'], 1)

    def test_positive_internal_cost_and_hard_forced_global_reversal(self):
        raw = star({'a': [1, 2, 3], 'b': [5, 6]}, leaf_orders={'a': [0, 2, 1], 'b': [0, 1]})
        raw['orientation_constraints'] = [dict(a='A:center', b='A:center', xor=0),
                                         dict(a='A:center', b='B:b', xor=0)]
        result = self.check_exhaustive(raw)
        self.assertEqual(result['bundle_lower'], 1)
        self.assertEqual(result['upper'], 1)

    def test_component_orders_and_isolated_chromosomes(self):
        raw = star({'a': [1, 3], 'b': [2, 4]})
        raw['species'][0]['chromosomes'].append(dict(id='isolate', length=20, blocks=[]))
        raw['species'][0]['chromosomes'].append(dict(id='extra', length=20,
            blocks=[dict(occurrence_id='extraA',homology_id='extra',start=1,end=2,strand='+')]))
        raw['species'][1]['chromosomes'].append(dict(id='extra', length=20,
            blocks=[dict(occurrence_id='extraB',homology_id='extra',start=1,end=2,strand='+')]))
        self.check_exhaustive(raw)

    def test_center_reversal_reverses_leaf_order_and_preserves_formula(self):
        raw = star({'a': [1, 2], 'b': [4, 5]}, reversed_center_row=True,
                   leaf_orders={'a': [1, 0], 'b': [0, 1]})
        result = self.check_exhaustive(raw)
        self.assertEqual(result['layout']['chromosome_orientation']['A:center'], -1)
        self.assertEqual(result['layout']['chromosome_order']['B'], ['b', 'a'])

    def test_opposite_direction_stars_in_same_interface(self):
        raw = star({'a': [1, 3], 'b': [2, 4]})
        other = star({'x': [1, 3], 'y': [2, 4]})
        for row in other['species']:
            destination = raw['species'][1 if row['id'] == 'A' else 0]
            for c in row['chromosomes']:
                c['id'] = 'second_' + c['id']
                destination['chromosomes'].append(c)
        self.check_exhaustive(raw)

    def test_all_three_partner_words_with_two_links_each(self):
        words = sorted(set(itertools.permutations('aabbcc')))
        methods = set()
        for word in words:
            with self.subTest(word=''.join(word)):
                positions = {v: [i+1 for i, x in enumerate(word) if x == v] for v in 'abc'}
                result = self.check_exhaustive(star(positions))
                methods.add(result['components'][0]['method'])
        self.assertIn('ACYCLIC_PAIR_PREFERENCES', methods)

    def test_subset_table_against_every_permutation(self):
        # Includes ties and a preference cycle; recurrence itself uses integers.
        costs = [[0, 5, 4, 2], [4, 0, 5, 2], [5, 4, 0, 1], [2, 2, 3, 0]]
        order, table = _subset_dp(costs)
        for mask in range(1 << 4):
            vertices = [v for v in range(4) if mask & (1 << v)]
            expected = min(sum(costs[a][b] for i,a in enumerate(p) for b in p[i+1:])
                           for p in itertools.permutations(vertices))
            self.assertEqual(table[mask], expected)
        self.assertEqual(sum(costs[a][b] for i,a in enumerate(order) for b in order[i+1:]), table[-1])

    def test_declines_nonstar_multilayer_and_conflicting_preferences(self):
        from test_structural import square
        self.assertEqual(construct(Evidence(square()))['status'], 'OUTSIDE_STAR_FOREST_CLASS')
        raw = star({'a': [1, 2]})
        raw['orientation_constraints'] = [dict(a='A:center', b='B:a', xor=1)]
        self.assertEqual(construct(Evidence(raw))['status'], 'OUTSIDE_COMPATIBLE_BUNDLE_PREFERENCES')
        raw = star({'a': [1, 2]})
        raw['species'].append(dict(id='C',chromosomes=[dict(id='1',length=20,blocks=[])]))
        self.assertEqual(construct(Evidence(raw))['status'], 'OUTSIDE_TWO_ROW_CLASS')


if __name__ == '__main__':
    unittest.main()
