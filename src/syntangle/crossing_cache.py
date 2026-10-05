"""Exact cached crossing terms for whole-chromosome local moves.

Each pair of homology links depends on at most two binary decisions: relative
chromosome order, or orientation when both anchors are on the same chromosome.
Terms aggregate identical decisions; no homology or internal order is changed.
"""
from itertools import combinations

from .layout import AmbiguousHomologyError, _anchor_key, _occurrences_by_species_homology


def _relation(a, b):
    return (a > b) - (a < b)


def _decision(a, b):
    ca, ba = a
    cb, bb = b
    if ca.ref == cb.ref:
        relations = tuple(_relation(_anchor_key(ca, ba, 0, sign),
                                    _anchor_key(cb, bb, 0, sign)) for sign in (1, -1))
        return (ca.ref,), relations
    refs = tuple(sorted((ca.ref, cb.ref)))
    relation = -1 if ca.ref == refs[0] else 1
    return refs, (relation, -relation)


class CrossingCostCache:
    def __init__(self, fixture):
        index = _occurrences_by_species_homology(fixture)
        tables = {}
        for left, right in zip(fixture.species_ids, fixture.species_ids[1:]):
            homologies = sorted(set(index.get(left, {})) & set(index.get(right, {})))
            links = []
            for homology in homologies:
                a, b = index[left][homology], index[right][homology]
                if len(a) != 1 or len(b) != 1:
                    raise AmbiguousHomologyError(f"Homology {homology!r} is not one-to-one between {left} and {right}")
                links.append((a[0], b[0]))
            # Precompute side comparisons once for each chromosome/anchor pair.
            for a, b in combinations(links, 2):
                dl, rl = _decision(a[0], b[0])
                dr, rr = _decision(a[1], b[1])
                costs = tables.setdefault((dl, dr), [0, 0, 0, 0])
                for x in (0, 1):
                    for y in (0, 1):
                        costs[2*x+y] += int(rl[x] * rr[y] < 0)
        self.constant = 0
        self.terms = []
        self.by_decision = {}
        for (a, b), costs in tables.items():
            if len(set(costs)) == 1:
                self.constant += costs[0]
                continue
            i = len(self.terms)
            self.terms.append((a, b, tuple(costs)))
            self.by_decision.setdefault(a, set()).add(i)
            self.by_decision.setdefault(b, set()).add(i)
        self.decisions_by_ref = {}
        for decision in self.by_decision:
            for ref in decision:
                self.decisions_by_ref.setdefault(ref, set()).add(decision)

    @staticmethod
    def _ranks(state):
        return {ref: rank for refs in state.chromosome_order.values() for rank, ref in enumerate(refs)}

    @staticmethod
    def _value(decision, state, ranks):
        if len(decision) == 1:
            return int(state.chromosome_orientation[decision[0]] == -1)
        return int(ranks[decision[0]] > ranks[decision[1]])

    def prepare(self, state):
        self.state = state
        self.ranks = self._ranks(state)
        self.values = {d: self._value(d, state, self.ranks) for d in self.by_decision}
        self.costs = [c[2*self.values[a]+self.values[b]] for a, b, c in self.terms]
        self.total = self.constant + sum(self.costs)
        return self.total

    def score_candidate(self, candidate):
        ranks = self._ranks(candidate)
        changed_refs = {ref for ref in ranks if ranks[ref] != self.ranks[ref]
                        or candidate.chromosome_orientation[ref] != self.state.chromosome_orientation[ref]}
        decisions = set()
        for ref in changed_refs:
            decisions.update(self.decisions_by_ref.get(ref, ()))
        changed = {}
        affected = set()
        for decision in decisions:
            value = self._value(decision, candidate, ranks)
            if value != self.values[decision]:
                changed[decision] = value
                affected.update(self.by_decision[decision])
        delta = 0
        for i in affected:
            a, b, costs = self.terms[i]
            delta += costs[2*changed.get(a, self.values[a])+changed.get(b, self.values[b])] - self.costs[i]
        return self.total + delta
