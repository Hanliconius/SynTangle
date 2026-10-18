"""Exact row-state recurrence for any number of adjacent-linked species rows.

No production imports or numerical solvers. Formal scope: SPECIES_AND_SIZE_THEOREMS.md.
"""
import argparse
import copy
import itertools
import json
import math
from pathlib import Path
from proof import Evidence, read_json


def row_states(evidence, species, local_constraints):
    names = evidence.rows[species]
    refs = [species + ':' + name for name in names]
    states = []
    for signs in itertools.product((1, -1), repeat=len(refs)):
        flips = dict(zip(refs, signs))
        if any(((flips[a] == -1) ^ (flips[b] == -1)) != p
               for a, b, p in local_constraints):
            continue
        for order in itertools.permutations(names):
            states.append(dict(order=list(order), flips=flips))
    return states


def interface_cost(evidence, index, left, right):
    species = list(evidence.rows)
    ranks = {species[index] + ':' + name: i for i, name in enumerate(left['order'])}
    ranks.update({species[index+1] + ':' + name: i for i, name in enumerate(right['order'])})
    flips = {**left['flips'], **right['flips']}
    def key(anchor):
        ref, x = anchor
        return ranks[ref], x if flips[ref] == 1 else evidence.lengths[ref] - x
    cost = 0
    for (a, b), (c, d) in itertools.combinations(evidence.links[index], 2):
        x, y = (key(a) > key(c)) - (key(a) < key(c)), (key(b) > key(d)) - (key(b) < key(d))
        cost += x * y < 0
    return cost


def construct(evidence, max_row_states=512, max_transitions=2_000_000):
    species = list(evidence.rows)
    row_of = {sp + ':' + name: i for i, (sp, names) in enumerate(evidence.rows.items()) for name in names}
    local = [[] for _ in species]
    adjacent = [[] for _ in species[1:]]
    for a, b, p in evidence.constraints:
        i, j = row_of[a], row_of[b]
        if abs(i-j) > 1:
            return dict(status='OUTSIDE_LOCAL_ORIENTATION_CONSTRAINTS')
        (local[i] if i == j else adjacent[min(i, j)]).append((a, b, p))
    capacities = [math.factorial(len(evidence.rows[sp])) * 2**len(evidence.rows[sp]) for sp in species]
    if max(capacities) > max_row_states:
        return dict(status='UNRESOLVED_ROW_STATE_LIMIT', unconstrained_row_states=capacities)
    states = [row_states(evidence, sp, eq) for sp, eq in zip(species, local)]
    sizes = list(map(len, states))
    transitions = sum(a*b for a, b in itertools.pairwise(sizes))
    if transitions > max_transitions:
        return dict(status='UNRESOLVED_TRANSITION_LIMIT', transitions=transitions)
    if any(not row for row in states):
        return dict(status='INFEASIBLE_LOCAL_EQUATIONS')
    tables, parents = [[0] * sizes[0]], []
    interface_minima = []
    for t in range(len(species)-1):
        costs = [[None] * sizes[t+1] for _ in states[t]]
        for i, left in enumerate(states[t]):
            for j, right in enumerate(states[t+1]):
                flips = {**left['flips'], **right['flips']}
                if all(((flips[a] == -1) ^ (flips[b] == -1)) == p for a, b, p in adjacent[t]):
                    costs[i][j] = interface_cost(evidence, t, left, right)
        finite = [c for row in costs for c in row if c is not None]
        interface_minima.append(min(finite) if finite else None)
        table, parent = [], []
        for j in range(sizes[t+1]):
            options = [(value + costs[i][j], i) for i, value in enumerate(tables[-1])
                       if value is not None and costs[i][j] is not None]
            value, prior = min(options) if options else (None, None)
            table.append(value)
            parent.append(prior)
        tables.append(table)
        parents.append(parent)
    options = [(v, i) for i, v in enumerate(tables[-1]) if v is not None]
    if not options:
        return dict(status='INFEASIBLE_LOCAL_EQUATIONS')
    optimum, index = min(options)
    path = [index]
    for parent in reversed(parents):
        index = parent[index]
        path.append(index)
    path.reverse()
    selected = [row[i] for row, i in zip(states, path)]
    layout = dict(chromosome_order={sp: s['order'] for sp, s in zip(species, selected)},
                  chromosome_orientation={r: sign for s in selected for r, sign in s['flips'].items()})
    assert evidence.score(layout) == optimum
    return dict(status='EXACT_SPECIES_CHAIN_RECURRENCE', upper=optimum, lower=optimum,
                row_state_counts=sizes, transitions=transitions,
                independent_interface_minima=interface_minima,
                interface_minimum_sum=sum(interface_minima),
                compatibility_penalty=optimum-sum(interface_minima),
                layout=layout, selected_state_indices=path, prefix_tables=tables)


def construct_with_separators(raw, max_row_states=512, max_transitions=2_000_000,
                              max_cyclic_size=12):
    """Glue exact contiguous subproblems at singleton rows using mirror symmetry."""
    from fusion_theorem import construct as fusion_construct
    evidence = Evidence(raw)
    species = list(evidence.rows)
    row_of = {r: i for i, sp in enumerate(species) for r in evidence.lengths
              if r.split(':', 1)[0] == sp}
    if any(abs(row_of[a]-row_of[b]) > 1 for a,b,_ in evidence.constraints):
        return dict(status='OUTSIDE_LOCAL_ORIENTATION_CONSTRAINTS')
    if len(species) == 1:
        return construct(evidence, max_row_states, max_transitions)
    cuts = [0] + [i for i in range(1, len(species)-1)
                  if len(evidence.rows[species[i]]) == 1] + [len(species)-1]
    orders, flips, pieces, total = {}, {}, [], 0
    for a, b in itertools.pairwise(cuts):
        part = copy.deepcopy(raw)
        part['species'] = part['species'][a:b+1]
        members = {sp + ':' + name for sp in species[a:b+1] for name in evidence.rows[sp]}
        part['orientation_constraints'] = [dict(a=x,b=y,xor=p) for x,y,p in evidence.constraints
                                           if x in members and y in members]
        local = Evidence(part)
        result = fusion_construct(local, max_cyclic_size) if b-a == 1 else {}
        if 'layout' not in result:
            result = construct(local, max_row_states, max_transitions)
        if 'layout' not in result:
            return dict(status='UNRESOLVED_SEPARATOR_PIECE', rows=[a,b], detail=result)
        state = copy.deepcopy(result['layout'])
        shared = members & flips.keys()
        mirrored = bool(shared and flips[next(iter(shared))] != state['chromosome_orientation'][next(iter(shared))])
        if shared:
            assert len(shared) == 1
        if mirrored:
            state['chromosome_order'] = {sp: list(reversed(order)) for sp,order in state['chromosome_order'].items()}
            state['chromosome_orientation'] = {r: -s for r,s in state['chromosome_orientation'].items()}
        assert local.score(state) == result['upper']
        for sp, order in state['chromosome_order'].items():
            if sp in orders:
                assert orders[sp] == order
            orders[sp] = order
        flips.update(state['chromosome_orientation'])
        total += result['upper']
        pieces.append(dict(rows=[a,b], mirrored=mirrored, result=result))
    layout = dict(chromosome_order=orders,chromosome_orientation=flips)
    assert evidence.score(layout) == total
    return dict(status='EXACT_SINGLETON_SEPARATOR_CONSTRUCTION',upper=total,lower=total,
                layout=layout,pieces=pieces,separator_rows=cuts[1:-1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fixture', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--max-row-states', type=int, default=512)
    parser.add_argument('--max-transitions', type=int, default=2_000_000)
    parser.add_argument('--separators', action='store_true', help='Split at singleton rows and use fusion-star pieces when possible')
    parser.add_argument('--max-cyclic-size', type=int, default=12)
    args = parser.parse_args()
    if min(args.max_row_states, args.max_transitions) < 1:
        parser.error('State and transition caps must be positive')
    if args.max_cyclic_size < 0:
        parser.error('--max-cyclic-size must be nonnegative')
    raw = read_json(args.fixture)
    result = (construct_with_separators(raw,args.max_row_states,args.max_transitions,args.max_cyclic_size)
              if args.separators else construct(Evidence(raw),args.max_row_states,args.max_transitions))
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], 'U=', result.get('upper'), 'L=', result.get('lower'),
          'compatibility_penalty=', result.get('compatibility_penalty'))


if __name__ == '__main__':
    main()
