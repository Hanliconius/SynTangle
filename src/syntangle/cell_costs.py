"""Lossless crossing-factor aggregation by chromosome-pair cells.

Counts every supplied link pair, including endpoint ties. No anchor sampling,
biological contraction, or hard ordering assumptions are made.
"""
from bisect import bisect_left, bisect_right
from collections import defaultdict
from itertools import combinations_with_replacement
import time

from .crossing_cache import _decision
from .incidence import chromosome_node_id
from .layout import _anchor_key, _strict_inversion_count, _occurrences_by_species_homology, AmbiguousHomologyError


class PreparationDeadline(TimeoutError):
    """Preparation stopped before a complete model/search frontier existed."""


def check_deadline(deadline):
    if deadline is not None and time.perf_counter() >= deadline:
        raise PreparationDeadline('Crossing-factor preparation deadline reached')


def relation_counts(first, second):
    """Number of (a,b) pairs with a<b and a>b; equal endpoints contribute zero."""
    ordered = sorted(second)
    less = sum(len(ordered)-bisect_right(ordered, x) for x in first)
    greater = sum(bisect_left(ordered, x) for x in first)
    return less, greater


def iter_cell_costs(fixture, nodes, deadline=None):
    index = _occurrences_by_species_homology(fixture)
    for left, right in zip(fixture.species_ids, fixture.species_ids[1:]):
        check_deadline(deadline)
        cells = defaultdict(list)
        for h in sorted(set(index.get(left, {})) & set(index.get(right, {}))):
            a, b = index[left][h], index[right][h]
            if len(a) != 1 or len(b) != 1:
                raise AmbiguousHomologyError('Cell factors require explicit one-to-one link identities')
            if all(chromosome_node_id(c.ref) in nodes for c, _ in (a[0], b[0])):
                cells[a[0][0].ref, b[0][0].ref].append((a[0], b[0]))
        for key_a, key_b in combinations_with_replacement(sorted(cells), 2):
            check_deadline(deadline)
            first, second = cells[key_a], cells[key_b]
            if key_a == key_b:
                dl, dr = (key_a[0],), (key_a[1],)
                costs = []
                for x in (1, -1):
                    for y in (1, -1):
                        costs.append(_strict_inversion_count([
                            (_anchor_key(ca, ba, 0, x), _anchor_key(cb, bb, 0, y))
                            for (ca, ba), (cb, bb) in first]))
            else:
                dl, rl = _decision(first[0][0], second[0][0])
                dr, rr = _decision(first[0][1], second[0][1])
                if key_a[0] != key_b[0] and key_a[1] != key_b[1]:
                    weight = len(first)*len(second)
                    costs = [weight*int(rl[x]*rr[y] < 0) for x in (0, 1) for y in (0, 1)]
                else:
                    side = 0 if key_a[0] == key_b[0] else 1
                    counts = []
                    for sign in (1, -1):
                        def values(links):
                            return [_anchor_key(c, b, 0, sign)[1] for link in links for c,b in (link[side],)]
                        counts.append(relation_counts(values(first), values(second)))
                    costs = []
                    for x in (0, 1):
                        for y in (0, 1):
                            less, greater = counts[x if side == 0 else y]
                            other_relation = rr[y] if side == 0 else rl[x]
                            costs.append(less if other_relation > 0 else greater)
            yield dl, dr, tuple(costs)
