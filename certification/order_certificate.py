"""Independent full-space branch certificates with exact rational LP-bound replay.

SciPy proposes multipliers only; verification uses the standard library and never
trusts an LP status, floating objective, or production optimizer bound.
"""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import time
from proof import Evidence, IntegerModel, read_json, digest
from structural import lower_bound as structural_bound


def canonical_rows(model):
    # All rows become <= inequalities or equalities, in deterministic order.
    result = []
    for terms, sense, rhs in model.rows:
        sign = -1 if sense == '>=' else 1
        result.append(({v: sign*c for v,c in terms.items()},
                       '=' if sense == '=' else '<=', sign*rhs))
    return result


def bound(model, rows, fixed, multipliers):
    """For Ax<=b, lambda>=0: C >= c0-lambda*b+min_box(c+A^T lambda)x.

Equality multipliers may have either sign. No stationarity assumption is needed;
exactly accounting for the residual makes rounded multipliers safe.
"""
    if len(multipliers) != len(rows):
        raise ValueError('Multiplier count differs from original model')
    if any(v not in model.binary or type(b) is not int or b not in (0,1)
           for v,b in fixed.items()):
        raise ValueError('Invalid fixed binary decision')
    residual = {v: Fraction(model.objective.get(v,0))
                for v in model.binary + list(model.products.values())}
    lower = Fraction(model.constant)
    for (terms,sense,rhs), value in zip(rows,multipliers):
        value = Fraction(value)
        if sense == '<=' and value < 0:
            raise ValueError('Negative inequality multiplier')
        lower -= value*rhs
        for v,c in terms.items():
            residual[v] += value*c
    for v,c in residual.items():
        lower += c*fixed[v] if v in fixed else min(0,c)
    return max(0, math.ceil(lower))


def bindings(fixture, layout):
    return dict(fixture_sha256=digest(fixture), layout_sha256=digest(layout),
                model_checker_sha256=digest(Path(__file__).with_name('proof.py')),
                branch_checker_sha256=digest(__file__),
                structural_checker_sha256=digest(Path(__file__).with_name('structural.py')))


def replay(model, upper, certificate):
    rows = canonical_rows(model)
    nodes = certificate['nodes']
    seen = set()
    def visit(index, expected):
        if type(index) is not int or not 0 <= index < len(nodes) or index in seen:
            raise ValueError('Missing, shared or cyclic branch node')
        seen.add(index)
        node = nodes[index]
        if node['fixed'] != expected:
            raise ValueError('Branch assignment is not inherited exactly')
        local = bound(model,rows,expected,node['multipliers'])
        if node.get('lower') != local:
            raise ValueError('Rational bound replay differs')
        if 'children' in node:
            variable = node['branch']
            if variable not in model.binary or variable in expected or len(node['children']) != 2:
                raise ValueError('Invalid branch partition')
            child_bounds = [visit(child,{**expected,variable:bit})
                            for bit,child in enumerate(node['children'])]
            return max(local,min(child_bounds))
        return local
    if not nodes:
        raise ValueError('Missing root')
    lower = visit(0,{})
    if len(seen) != len(nodes):
        raise ValueError('Unreachable certificate nodes')
    if lower > upper:
        raise ValueError('Bound exceeds independently validated candidate')
    return dict(upper=upper,lower=lower,gap=upper-lower,
                status='EXACT_ORDER_CERTIFICATE' if lower == upper else 'UNRESOLVED_ORDER_GAP',
                nodes=len(nodes))


def generate(model, upper, seconds, max_nodes):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    rows = canonical_rows(model)
    variables = model.binary + list(model.products.values())
    indices = {v:i for i,v in enumerate(variables)}
    groups = [[i for i,(_,sense,_) in enumerate(rows) if sense == target]
              for target in ('<=','=')]
    def matrix(group):
        triples = [(j,indices[v],c) for j,i in enumerate(group) for v,c in rows[i][0].items()]
        if not triples:
            return coo_matrix((len(group),len(variables))).tocsr()
        a,b,c = zip(*triples)
        return coo_matrix((c,(a,b)),shape=(len(group),len(variables))).tocsr()
    matrices = [matrix(g) for g in groups]
    rhs = [np.array([rows[i][2] for i in g],dtype=float) for g in groups]
    costs = np.array([model.objective.get(v,0) for v in variables],dtype=float)
    deadline = time.monotonic()+seconds
    nodes = []
    pending = []
    def append(fixed):
        index = len(nodes)
        nodes.append(dict(fixed=fixed,multipliers=['0']*len(rows),
                          lower=bound(model,rows,fixed,['0']*len(rows))))
        pending.append(index)
        return index
    append({})
    solved = 0
    while pending and time.monotonic()<deadline:
        index = pending.pop()
        node = nodes[index]
        limits = [(node['fixed'].get(v,0),node['fixed'].get(v,1)) for v in variables]
        remaining = deadline-time.monotonic()
        if remaining <= 0:
            break
        result = linprog(costs,A_ub=matrices[0] if groups[0] else None,
                         b_ub=rhs[0] if groups[0] else None,
                         A_eq=matrices[1] if groups[1] else None,
                         b_eq=rhs[1] if groups[1] else None,bounds=limits,
                         method='highs',options={'time_limit':min(10.0,remaining),'threads':1})
        solved += 1
        # Infeasible floating solver status is NOT a proof. Keep safe zero multipliers.
        if result.success:
            proposed = [Fraction(0)]*len(rows)
            for group,part in zip(groups,(result.ineqlin,result.eqlin)):
                for i,value in zip(group,part.marginals):
                    proposed[i] = Fraction(float(-value)).limit_denominator(1000000)
            lower = bound(model,rows,node['fixed'],proposed)
            if lower >= node['lower']:
                node['multipliers'] = [str(x) for x in proposed]
                node['lower'] = lower
        if node['lower'] >= upper:
            continue
        free = [v for v in model.binary if v not in node['fixed']]
        if not free or len(nodes)+2 > max_nodes or time.monotonic()>=deadline:
            continue
        # Branches partition the original full binary space. No heuristic exclusion.
        if result.success:
            variable = min(free,key=lambda v:(abs(result.x[indices[v]]-0.5),v))
        else:
            variable = free[0]
        node['branch'] = variable
        node['children'] = [append({**node['fixed'],variable:bit}) for bit in (0,1)]
    return dict(nodes=nodes,lp_calls=solved,deadline=time.monotonic()>=deadline,
                variables=len(variables),binary_variables=len(model.binary),rows=len(rows))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['check','verify'])
    p.add_argument('fixture',type=Path)
    p.add_argument('layout',type=Path)
    p.add_argument('certificate',type=Path)
    p.add_argument('--seconds',type=float,default=600)
    p.add_argument('--max-nodes',type=int,default=1001)
    a = p.parse_args()
    if a.seconds < 0 or a.max_nodes < 1:
        raise ValueError('Invalid budget')
    evidence = Evidence(read_json(a.fixture))
    candidate = read_json(a.layout)
    upper = evidence.score(candidate)
    model = IntegerModel(evidence)
    if model.evaluate(model.assignment(candidate)) != upper:
        raise ValueError('Original geometry/model disagreement')
    hashes = bindings(a.fixture,a.layout)
    base_lower = structural_bound(evidence)['lower']
    def combined(certificate):
        result = replay(model,upper,certificate)
        result['ordering_lower'] = result['lower']
        result['structural_lower'] = base_lower
        result['lower'] = max(base_lower,result['lower'])
        if result['lower'] > upper:
            raise ValueError('Combined bound exceeds candidate')
        result['gap'] = upper-result['lower']
        result['status'] = 'EXACT_ORDER_CERTIFICATE' if result['gap']==0 else 'UNRESOLVED_ORDER_GAP'
        return result
    if a.action == 'check':
        certificate = {**hashes,**generate(model,upper,a.seconds,a.max_nodes)}
        summary = combined(certificate)
        certificate['summary'] = summary
        a.certificate.write_text(json.dumps(certificate,separators=(',',':'))+'\n')
    else:
        certificate = json.loads(a.certificate.read_text())
        if any(certificate.get(k)!=v for k,v in hashes.items()):
            raise ValueError('Input or checker hash mismatch')
        summary = combined(certificate)
        if certificate['summary'] != summary:
            raise ValueError('Saved conclusion differs from replay')
    print(('VERIFIED_' if a.action=='verify' else '')+summary['status'],
          'U=',upper,'L=',summary['lower'],'gap=',summary['gap'],'nodes=',summary['nodes'])

if __name__ == '__main__':
    main()
