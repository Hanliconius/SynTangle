import itertools
import random
import unittest
from proof import Evidence
from fusion_theorem import construct, _ordered_sccs, _subset_dp
from test_fusion_theorem import star


class FusionScalingTests(unittest.TestCase):
    def test_scc_decomposition_matches_all_orders(self):
        rng = random.Random(104)
        for n in range(2, 7):
            for repeat in range(12):
                costs = [[0 if i==j else rng.randrange(6) for j in range(n)] for i in range(n)]
                parts = _ordered_sccs(costs)
                order = []
                for part in parts:
                    local, _ = _subset_dp([[costs[i][j] for j in part] for i in part])
                    order.extend(part[i] for i in local)
                score = sum(costs[a][b] for i,a in enumerate(order) for b in order[i+1:])
                best = min(sum(costs[a][b] for i,a in enumerate(p) for b in p[i+1:])
                           for p in itertools.permutations(range(n)))
                self.assertEqual(score, best)

    def test_many_partners_only_three_are_cyclic(self):
        # Extend the exact 3-partner cycle by 97 ordered, noninterleaving partners.
        positions = {'a':[2,4,9], 'b':[1,6,8], 'c':[3,5,7]}
        positions.update({f'z{i:03}':[10+i] for i in range(97)})
        raw = star(positions)
        for row in raw['species']:
            for chromosome in row['chromosomes']:
                chromosome['length'] = 1000
        result = construct(Evidence(raw), max_degree=3)
        self.assertEqual(result['upper'], 13)
        item = result['components'][0]
        self.assertEqual(len(item['leaves']), 100)
        self.assertEqual(max(len(p['vertices']) for p in item['ordering_parts']), 3)
        self.assertEqual(sum(len(p['subset_table']) for p in item['ordering_parts']), 8+97*2)


if __name__ == '__main__':
    unittest.main()
