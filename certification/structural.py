"""Structural lower-bound certificates; independent of optimizer and MIP solver."""
import argparse
from collections import defaultdict, deque
import itertools as it
import json
from pathlib import Path
from proof import Evidence, read_json, digest


def bundles(evidence):
    result = []
    for links in evidence.links:
        groups = defaultdict(list)
        for left, right in links:
            groups[left[0], right[0]].append((left[1], right[1]))
        for (a,b), anchors in sorted(groups.items()):
            concordant = discordant = 0
            for (x,y), (u,v) in it.combinations(anchors,2):
                product = (x-u)*(y-v)
                concordant += product > 0
                discordant += product < 0
            result.append(dict(a=a,b=b,concordant=concordant,discordant=discordant))
    return result


def path(forest, start, end):
    """Unique forest path, including the empty path from a vertex to itself."""
    queue, previous = deque([start]), {start:None}
    while queue:
        node = queue.popleft()
        if node == end:
            edges = []
            while previous[node] is not None:
                node, edge = previous[node]
                edges.append(edge)
            return edges
        for neighbor, edge in forest[node]:
            if neighbor not in previous:
                previous[neighbor] = (node,edge)
                queue.append(neighbor)
    return None


def lower_bound(evidence):
    """B0 plus soft-edge-disjoint inconsistent parity cycle witnesses.

Hard edges may be reused between witnesses: they never pay a penalty. Soft
edges may not be reused, preventing double counting. Greedy packing is not
claimed to find a strongest bound or minimum frustration.
"""
    counts = bundles(evidence)
    base = sum(min(g['concordant'],g['discordant']) for g in counts)
    forest = defaultdict(list)
    edges, cycles, used = [], [], set()
    hard = [dict(a=a,b=b,parity=p,weight=None) for a,b,p in evidence.constraints]
    soft = [dict(a=g['a'],b=g['b'],parity=int(g['discordant']>g['concordant']),
                 weight=abs(g['discordant']-g['concordant']),bundle=i)
            for i,g in enumerate(counts) if g['discordant'] != g['concordant']]
    # Hard forest first. Equal-cost bundles impose no preferred parity.
    for edge in hard + sorted(soft,key=lambda e:(-e['weight'],e['bundle'])):
        index = len(edges)
        edges.append(edge)
        found = path(forest,edge['a'],edge['b'])
        if found is None:
            forest[edge['a']].append((edge['b'],index))
            forest[edge['b']].append((edge['a'],index))
            continue
        parity = edge['parity']
        for i in found:
            parity ^= edges[i]['parity']
        if not parity:
            continue
        witness = found + [index]
        paid = {i for i in witness if edges[i]['weight'] is not None}
        if not paid:
            raise ValueError('Inconsistent hard orientation equations: no legal layout')
        if paid & used:
            continue
        penalty = min(edges[i]['weight'] for i in paid)
        cycles.append(dict(edges=witness,penalty=penalty))
        used.update(paid)
    bound = base + sum(c['penalty'] for c in cycles)
    return dict(base=base,lower=bound,bundles=counts,parity_edges=edges,cycles=cycles)


def construct_thin(evidence):
    """Construct an optimum for row-thin incidence components if preferences fit.

Uses all scored links, not merely nonzero orientation-cost edges. No evidence
is removed. Declines unsupported structures instead of claiming impossibility.
"""
    incidence = defaultdict(set)
    for links in evidence.links:
        for (a,_),(b,_) in links:
            incidence[a].add(b)
            incidence[b].add(a)
    components, assigned = [], set()
    for ref in evidence.lengths:
        if ref in assigned:
            continue
        members, queue = set(), [ref]
        while queue:
            node = queue.pop()
            if node in members:
                continue
            members.add(node)
            queue.extend(incidence[node]-members)
        assigned.update(members)
        species = [r.split(':',1)[0] for r in members]
        if len(species) != len(set(species)):
            return dict(status='OUTSIDE_THIN_CLASS')
        components.append(members)
    equations = defaultdict(list)
    edges = [(a,b,p) for a,b,p in evidence.constraints]
    for g in bundles(evidence):
        if g['concordant'] != g['discordant']:
            edges.append((g['a'],g['b'],int(g['discordant']>g['concordant'])))
    for a,b,p in edges:
        equations[a].append((b,p))
        equations[b].append((a,p))
    bits = {}
    for ref in evidence.lengths:
        if ref in bits:
            continue
        bits[ref] = 0
        queue = [ref]
        while queue:
            node = queue.pop()
            for other,p in equations[node]:
                value = bits[node]^p
                if other in bits:
                    if bits[other] != value:
                        return dict(status='PREFERENCES_INCONSISTENT')
                else:
                    bits[other] = value
                    queue.append(other)
    components.sort(key=lambda members:min(members))
    index = {r:i for i,members in enumerate(components) for r in members}
    candidate = dict(chromosome_order={sp:sorted(names,key=lambda n:index[sp+':'+n])
                                     for sp,names in evidence.rows.items()},
                     chromosome_orientation={r:1-2*b for r,b in bits.items()})
    score = evidence.score(candidate)
    bound = lower_bound(evidence)['lower']
    if score != bound:
        raise ValueError('Construction/theorem disagreement')
    return dict(status='EXACT_THIN_CONSTRUCTION',upper=score,lower=bound,layout=candidate)


def report(evidence,candidate):
    bound = lower_bound(evidence)
    score = evidence.score(candidate)
    if score < bound['lower']:
        raise ValueError('Structural lower bound exceeds independent score')
    return {**bound,'upper':score,'gap':score-bound['lower'],
            'status':'EXACT_STRUCTURAL' if score == bound['lower'] else 'UNRESOLVED_STRUCTURAL_GAP'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['check','verify','construct'])
    parser.add_argument('fixture',type=Path)
    parser.add_argument('layout',type=Path,help='Candidate JSON; construction output for construct')
    parser.add_argument('certificate',type=Path)
    args = parser.parse_args()
    evidence = Evidence(read_json(args.fixture))
    if args.action == 'construct':
        constructed = construct_thin(evidence)
        if 'layout' not in constructed:
            print(constructed['status'])
            return
        args.layout.write_text(json.dumps(constructed['layout'],indent=2)+'\n')
    certificate = {**report(evidence,read_json(args.layout)),
        'fixture_sha256':digest(args.fixture),'layout_sha256':digest(args.layout),
        'structural_checker_sha256':digest(__file__),
        'geometry_checker_sha256':digest(Path(__file__).with_name('proof.py'))}
    if args.action == 'verify':
        saved = read_json(args.certificate)
        if saved != certificate or certificate['status'] != 'EXACT_STRUCTURAL':
            raise ValueError('Replay did not validate an exact structural certificate')
        print('VERIFIED_EXACT_STRUCTURAL',certificate['upper'])
    else:
        args.certificate.write_text(json.dumps(certificate,indent=2)+'\n')
        print(certificate['status'],'U=',certificate['upper'],'L=',certificate['lower'])


if __name__ == '__main__':
    main()
