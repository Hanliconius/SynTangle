"""Independent exact certification experiment. Standard library only.

No imports from syntangle: original JSON decimal coordinates are rationals.
The exhaustive certificate is verified by replay, not by trusting a result flag.
The LP export is NOT itself a certificate.
"""
import argparse
import hashlib
import itertools as it
import json
import math
import time
from fractions import Fraction
from pathlib import Path


def read_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError('Duplicate JSON key: ' + key)
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), parse_float=Fraction,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)),
                      object_pairs_hook=unique)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Evidence:
    def __init__(self, raw):
        if raw.get('fixture_version') != 1:
            raise ValueError('Only fixture version 1 is supported')
        supported = {'fixture_version', 'id', 'title', 'purpose', 'species',
                     'orientation_constraints', 'expected', 'metadata', 'provenance'}
        if set(raw) - supported:
            raise ValueError('Unsupported top-level fields: ' + str(set(raw)-supported))
        self.rows, self.lengths, self.occurrences = {}, {}, {}
        seen = set()
        for row in raw['species']:
            sp = row['id']
            if not isinstance(sp, str) or ':' in sp or sp in self.rows:
                raise ValueError('Invalid or duplicate species id')
            self.rows[sp] = []
            self.occurrences[sp] = {}
            for chrom in row['chromosomes']:
                if set(chrom) - {'id', 'length', 'display_rank', 'display_orientation', 'blocks'}:
                    raise ValueError('Unsupported chromosome fields')
                name = chrom['id']
                if not isinstance(name, str) or name in self.rows[sp]:
                    raise ValueError('Invalid or duplicate chromosome id')
                ref = sp + ':' + name
                length = Fraction(chrom['length'])
                if length <= 0:
                    raise ValueError('Invalid chromosome length')
                self.rows[sp].append(name)
                self.lengths[ref] = length
                for block in chrom['blocks']:
                    if set(block) - {'occurrence_id', 'homology_id', 'start', 'end', 'strand'}:
                        raise ValueError('Unsupported block fields; weighted objectives unsupported')
                    occurrence, h = block['occurrence_id'], block['homology_id']
                    if occurrence in seen or h in self.occurrences[sp]:
                        raise ValueError('Duplicate occurrence or ambiguous homology')
                    seen.add(occurrence)
                    start, end = Fraction(block['start']), Fraction(block['end'])
                    if not 0 <= start < end <= length:
                        raise ValueError('Invalid block interval')
                    self.occurrences[sp][h] = (ref, (start + end) / 2)
        if not self.rows or any(not row for row in self.rows.values()):
            raise ValueError('Empty species rows unsupported')
        self.constraints = []
        for c in raw.get('orientation_constraints', []):
            if set(c) != {'a', 'b', 'xor'} or c['a'] not in self.lengths or c['b'] not in self.lengths or type(c['xor']) is not int or c['xor'] not in (0, 1):
                raise ValueError('Invalid orientation equation')
            self.constraints.append((c['a'], c['b'], c['xor']))
        self.links = []
        for left, right in it.pairwise(self.rows):
            common = sorted(self.occurrences[left].keys() & self.occurrences[right].keys())
            self.links.append([(self.occurrences[left][h], self.occurrences[right][h]) for h in common])

    def validate(self, state):
        orders, flips = state['chromosome_order'], state['chromosome_orientation']
        if set(orders) != set(self.rows) or set(flips) != set(self.lengths):
            raise ValueError('Layout must cover all species and chromosomes')
        for sp, names in self.rows.items():
            if len(orders[sp]) != len(names) or set(orders[sp]) != set(names):
                raise ValueError('Layout is not a whole-chromosome permutation')
        if any(type(s) is not int or s not in (-1, 1) for s in flips.values()):
            raise ValueError('Orientation must be integer +1 or -1')
        for a, b, parity in self.constraints:
            if ((flips[a] == -1) ^ (flips[b] == -1)) != parity:
                raise ValueError('Hard orientation constraint violated')

    def score(self, state):
        self.validate(state)
        ranks = {sp+':'+name: i for sp, names in state['chromosome_order'].items() for i, name in enumerate(names)}
        def key(anchor):
            ref, midpoint = anchor
            coordinate = midpoint if state['chromosome_orientation'][ref] == 1 else self.lengths[ref] - midpoint
            return ranks[ref], coordinate
        total = 0
        for links in self.links:
            for (a, b), (c, d) in it.combinations(links, 2):
                # Strict inequality: shared/tied anchors do not cross.
                left, right = (key(a) > key(c)) - (key(a) < key(c)), (key(b) > key(d)) - (key(b) < key(d))
                total += left * right < 0
        return total

    def space_size(self):
        return math.prod(math.factorial(len(r)) for r in self.rows.values()) * 2**len(self.lengths)

    def states(self):
        refs = list(self.lengths)
        def orders(index, current):
            if index == len(self.rows):
                yield dict(current)
                return
            sp = list(self.rows)[index]
            for perm in it.permutations(self.rows[sp]):
                current[sp] = list(perm)
                yield from orders(index+1, current)
        for bits in it.product((1, -1), repeat=len(refs)):
            flips = dict(zip(refs, bits))
            if any(((flips[a] == -1) ^ (flips[b] == -1)) != p for a, b, p in self.constraints):
                continue
            for order in orders(0, {}):
                yield {'chromosome_order': order, 'chromosome_orientation': flips}


def exhaustive(evidence, candidate, max_states, seconds):
    upper = evidence.score(candidate)
    if upper == 0:
        return dict(status='EXACT_ZERO', upper=0, lower=0, checked=0)
    if evidence.space_size() > max_states:
        return dict(status='UNRESOLVED_SPACE_LIMIT', upper=upper, lower=0, checked=0)
    deadline, checked = time.monotonic() + seconds, 0
    for state in evidence.states():
        if time.monotonic() >= deadline:
            return dict(status='UNRESOLVED_DEADLINE', upper=upper, lower=0, checked=checked)
        score = evidence.score(state)
        checked += 1
        if score < upper:
            return dict(status='COUNTEREXAMPLE', upper=upper, lower=0, checked=checked,
                        better_score=score, better_layout=state)
    return dict(status='EXACT_EXHAUSTIVE', upper=upper, lower=upper, checked=checked)


class IntegerModel:
    """Full-space model: no decomposition, symmetry fixing, or imported reductions.

Each crossing is the XOR of two endpoint precedence predicates. Products are
deduplicated and linearized exactly. Variables x_ab=1 mean a occurs after b.
"""
    def __init__(self, evidence):
        self.binary, self.rows, self.products = [], [], {}
        self.constant, self.objective = 0, {}
        self.flips = {ref: self.var() for ref in evidence.lengths}
        self.orders = {}
        for sp, names in evidence.rows.items():
            refs = [sp+':'+n for n in names]
            for a, b in it.combinations(refs, 2):
                self.orders[a, b] = self.var()
            for a, b, c in it.combinations(refs, 3):
                terms = {self.orders[a,b]: 1, self.orders[b,c]: 1, self.orders[a,c]: -1}
                self.rows += [(terms, '>=', 0), (terms, '<=', 1)]
        for a, b, p in evidence.constraints:
            u, v = self.flips[a], self.flips[b]
            if u == v:
                if p:
                    raise ValueError('Contradictory self orientation constraint')
            else:
                self.rows.append(({u: 1, v: 1 if p else -1}, '=', p))
        def predicate(a, b):
            r, m = a
            s, n = b
            if r == s:
                if m == n:
                    return None
                return (self.flips[r], int(m > n), -1 if m > n else 1)
            if (r,s) in self.orders:
                return self.orders[r,s], 0, 1
            return self.orders[s,r], 1, -1
        for links in evidence.links:
            for (a,b), (c,d) in it.combinations(links, 2):
                left, right = predicate(a,c), predicate(b,d)
                if left is None or right is None:
                    continue
                variables = sorted({left[0], right[0]})
                def cost(bits):
                    return int((left[1]+left[2]*bits[left[0]]) != (right[1]+right[2]*bits[right[0]]))
                zero = dict.fromkeys(variables, 0)
                base = cost(zero)
                self.constant += base
                increments = []
                for v in variables:
                    delta = cost({**zero, v:1}) - base
                    self.add(v, delta)
                    increments.append(delta)
                if len(variables) == 2:
                    q = cost(dict.fromkeys(variables, 1)) - base - sum(increments)
                    if q:
                        self.add(self.product(*variables), q)

    def var(self):
        name = 'x' + str(len(self.binary))
        self.binary.append(name)
        return name

    def add(self, variable, coefficient):
        self.objective[variable] = self.objective.get(variable, 0) + coefficient

    def product(self, a, b):
        pair = tuple(sorted((a,b)))
        if pair not in self.products:
            z = 'z' + str(len(self.products))
            self.products[pair] = z
            self.rows += [({z:1, a:-1}, '<=', 0), ({z:1, b:-1}, '<=', 0),
                          ({z:1, a:-1, b:-1}, '>=', -1)]
        return self.products[pair]

    def assignment(self, state):
        ranks = {sp+':'+n: i for sp, row in state['chromosome_order'].items() for i,n in enumerate(row)}
        values = {v:int(state['chromosome_orientation'][r] == -1) for r,v in self.flips.items()}
        values.update({v:int(ranks[a] > ranks[b]) for (a,b),v in self.orders.items()})
        values.update({z:values[a]*values[b] for (a,b),z in self.products.items()})
        return values

    def evaluate(self, values):
        return self.constant + sum(c*values[v] for v,c in self.objective.items())

    def lp(self, upper):
        def expression(terms):
            return ' '.join(('+' if c >= 0 else '-') + f' {abs(c)} {v}' for v,c in terms.items() if c) or '0 dummy'
        rows = self.rows + [(self.objective, '<=', upper - 1 - self.constant)]
        lines = ['Minimize', ' obj: dummy', 'Subject To']
        lines += [f' c{i}: {expression(terms)} {sense} {rhs}' for i,(terms,sense,rhs) in enumerate(rows)]
        lines += ['Bounds', ' dummy = 0'] + [f' 0 <= {z} <= 1' for z in self.products.values()]
        lines += ['Binary'] + [' '+v for v in self.binary] + ['End']
        return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['check', 'verify', 'export'])
    parser.add_argument('fixture', type=Path)
    parser.add_argument('layout', type=Path, help='Plain LayoutState.to_dict() JSON')
    parser.add_argument('output', type=Path, help='Certificate JSON or new export directory')
    parser.add_argument('--max-states', type=int, default=1000000)
    parser.add_argument('--seconds', type=float, default=600)
    args = parser.parse_args()
    evidence, candidate = Evidence(read_json(args.fixture)), read_json(args.layout)
    upper = evidence.score(candidate)
    hashes = {'fixture_sha256':digest(args.fixture), 'layout_sha256':digest(args.layout),
              'checker_sha256':digest(__file__)}
    if args.action == 'export':
        args.output.mkdir(parents=True, exist_ok=False)
        model = IntegerModel(evidence)
        if model.evaluate(model.assignment(candidate)) != upper:
            raise ValueError('Independent geometry/model disagreement')
        lp_path = args.output / 'no_better.lp'
        lp_path.write_text(model.lp(upper))
        (args.output/'manifest.json').write_text(json.dumps({**hashes, 'upper':upper,
            'model_sha256':digest(lp_path), 'status':'EXPORTED_NOT_CERTIFIED'}, indent=2)+'\n')
        # Paths relative to export directory; enable exact mode BEFORE read.
        (args.output/'scip.commands').write_text('set exact enable TRUE\nset certificate filename proof.vipr\nread no_better.lp\noptimize\nquit\n')
        print('EXPORTED_NOT_CERTIFIED', args.output)
        return
    if args.action == 'verify':
        saved = read_json(args.output)
        if any(saved.get(k) != v for k,v in hashes.items()):
            raise ValueError('Certificate input/checker fingerprint mismatch')
    result = exhaustive(evidence, candidate, args.max_states, args.seconds)
    if args.action == 'verify':
        if result['status'] not in ('EXACT_ZERO','EXACT_EXHAUSTIVE') or saved != {**hashes, **result}:
            raise ValueError('Certificate replay did not establish the claimed optimum')
        print('VERIFIED_BY_EXACT_REPLAY', result['upper'])
    else:
        args.output.write_text(json.dumps({**hashes, **result}, indent=2)+'\n')
        print(result['status'], 'U=', upper, 'checked=', result['checked'])


if __name__ == '__main__':
    main()
