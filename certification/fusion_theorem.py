"""Exact two-row fusion/fission star families; standard library, no optimizer.

The mathematical argument and scope are in FUSION_FISSION_THEOREMS.md.
OUTSIDE/UNRESOLVED statuses never assert non-optimality or infeasibility.
"""
import argparse
from collections import defaultdict
import itertools
import heapq
import json
from pathlib import Path

from proof import Evidence, read_json
from structural import bundles


def _components(evidence):
    neighbors = {r: set() for r in evidence.lengths}
    for links in evidence.links:
        for (a, _), (b, _) in links:
            neighbors[a].add(b)
            neighbors[b].add(a)
    groups, seen = [], set()
    for root in sorted(neighbors):
        if root in seen:
            continue
        members, todo = set(), [root]
        while todo:
            r = todo.pop()
            if r not in members:
                members.add(r)
                todo.extend(neighbors[r] - members)
        seen.update(members)
        center = min(members, key=lambda r: (-len(neighbors[r]), r))
        if any(r != center and len(neighbors[r]) != 1 for r in members):
            return None
        groups.append((center, sorted(members - {center})))
    return groups


def _preferred_flips(evidence, counts):
    equations = defaultdict(list)
    edges = list(evidence.constraints)
    for g in counts:
        if g['concordant'] != g['discordant']:
            edges.append((g['a'], g['b'], int(g['discordant'] > g['concordant'])))
    for a, b, parity in edges:
        equations[a].append((b, parity))
        equations[b].append((a, parity))
    bits = {}
    for root in evidence.lengths:
        if root in bits:
            continue
        bits[root] = 0
        todo = [root]
        while todo:
            r = todo.pop()
            for s, parity in equations[r]:
                bit = bits[r] ^ parity
                if s in bits and bits[s] != bit:
                    return None
                if s not in bits:
                    bits[s] = bit
                    todo.append(s)
    return {r: 1 - 2 * bit for r, bit in bits.items()}


def _topological(costs):
    n = len(costs)
    edges = [(i, j) for i in range(n) for j in range(n)
             if costs[i][j] < costs[j][i]]
    children, degree = defaultdict(list), [0] * n
    for a, b in edges:
        children[a].append(b)
        degree[b] += 1
    ready = [v for v in range(n) if degree[v] == 0]
    heapq.heapify(ready)
    order = []
    while ready:
        a = heapq.heappop(ready)
        order.append(a)
        for b in children[a]:
            degree[b] -= 1
            if degree[b] == 0:
                heapq.heappush(ready, b)
    return order if len(order) == n else None


def _subset_dp(costs):
    """D[S]: append v last, paying every preceding u -> v cost."""
    n = len(costs)
    table, last = [0] * (1 << n), [-1] * (1 << n)
    for mask in range(1, 1 << n):
        choices = []
        for v in range(n):
            if mask & (1 << v):
                prior = mask ^ (1 << v)
                choices.append((table[prior] + sum(costs[u][v] for u in range(n)
                                                   if prior & (1 << u)), v))
        table[mask], last[mask] = min(choices)
    order, mask = [], (1 << n) - 1
    while mask:
        v = last[mask]
        order.append(v)
        mask ^= 1 << v
    return order[::-1], table


def construct(evidence, max_degree=12):
    if len(evidence.rows) != 2:
        return {'status': 'OUTSIDE_TWO_ROW_CLASS'}
    groups = _components(evidence)
    if groups is None:
        return {'status': 'OUTSIDE_STAR_FOREST_CLASS'}
    counts = bundles(evidence)
    flips = _preferred_flips(evidence, counts)
    if flips is None:
        return {'status': 'OUTSIDE_COMPATIBLE_BUNDLE_PREFERENCES'}
    order = {sp: [] for sp in evidence.rows}
    base = sum(min(g['concordant'], g['discordant']) for g in counts)
    lower, witnesses = base, []
    for center, leaves in groups:
        anchors = {leaf: [] for leaf in leaves}
        for (a, x), (b, y) in evidence.links[0]:
            if a == center:
                anchors[b].append(x)
            elif b == center:
                anchors[a].append(y)
        costs = [[sum(x > y for x in anchors[a] for y in anchors[b])
                  if a != b else 0 for b in leaves] for a in leaves]
        pair_lower = sum(min(costs[i][j], costs[j][i])
                         for i, j in itertools.combinations(range(len(leaves)), 2))
        permutation = _topological(costs)
        table = None
        if permutation is not None:
            value, method = pair_lower, 'ACYCLIC_PAIR_PREFERENCES'
        elif len(leaves) <= max_degree:
            permutation, table = _subset_dp(costs)
            value, method = table[-1], 'EXACT_SUBSET_RECURRENCE'
        else:
            return {'status': 'UNRESOLVED_CYCLIC_DEGREE_LIMIT',
                    'center': center, 'degree': len(leaves), 'limit': max_degree}
        lower += value
        native_order = [leaves[i] for i in permutation]
        actual_order = native_order if flips[center] == 1 else native_order[::-1]
        for r in [center] + actual_order:
            sp, name = r.split(':', 1)
            order[sp].append(name)
        witnesses.append(dict(center=center, leaves=leaves, costs=costs,
                              native_leaf_order=native_order, method=method,
                              pair_lower=pair_lower, exact_between_bundle_cost=value,
                              subset_table=table))
    layout = dict(chromosome_order=order, chromosome_orientation=flips)
    upper = evidence.score(layout)
    if upper != lower:
        raise AssertionError('Construction does not attain the proved formula')
    return dict(status='EXACT_FUSION_STAR_CONSTRUCTION', upper=upper, lower=lower,
                bundle_lower=base, layout=layout, components=witnesses)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--max-degree', type=int, default=12,
                        help='Subset DP cap for cyclic preferences; acyclic stars have no cap')
    args = parser.parse_args()
    if args.max_degree < 0:
        parser.error('--max-degree must be nonnegative')
    result = construct(Evidence(read_json(args.fixture)), args.max_degree)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], 'U=', result.get('upper'), 'L=', result.get('lower'))


if __name__ == '__main__':
    main()
