"""Graph-selected exact four-cycle cuts and weakest-bound-first proof search.

Floating LPs suggest multipliers; the independent rational checker remains the
only authority. Optional cuts are valid on every binary assignment, not heuristic
restrictions on chromosome orders.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import heapq
import itertools
import json
from pathlib import Path
import time
from proof import Evidence,IntegerModel,read_json,digest
from structural import lower_bound as structural_bound
from order_certificate import bound,canonical_rows,replay


def cycle_cuts(model, limit, max_examined=100000):
    """Even parity around four-cycles of interacting binary decisions.

    d_uv = x_u + x_v - 2 z_uv is XOR on binary assignments. For every
    odd subset F of cycle edges, sum_F d - sum_not_F d <= |F|-1.
    A four-cycle has eight such inequalities. Every real assignment has
    even XOR parity, so every inequality is valid without excluding a layout.
    """
    if type(limit) is not int or limit<0:raise ValueError('Invalid cycle limit')
    if limit==0:return [],0
    products={tuple(sorted(pair)):v for pair,v in model.products.items()}
    neighbors=defaultdict(set)
    for a,b in products:neighbors[a].add(b);neighbors[b].add(a)
    opposites=defaultdict(set)
    for center in sorted(neighbors):
        for a,c in itertools.combinations(sorted(neighbors[center]),2):opposites[a,c].add(center)
    selected=[];examined=0
    for (a,c),middle in sorted(opposites.items()):
        for b,d in itertools.combinations(sorted(middle),2):
            if a!=min(a,b,c,d):continue  # each undirected cycle exactly once
            cycle=(a,b,c,d)
            edges=[tuple(sorted(edge)) for edge in ((a,b),(b,c),(c,d),(d,a))]
            weight=sum(abs(model.objective.get(products[e],0)) for e in edges)
            item=(weight,cycle)
            if len(selected)<limit:heapq.heappush(selected,item)
            elif item>selected[0]:heapq.heapreplace(selected,item)
            examined+=1
            if examined>=max_examined:break
        if examined>=max_examined:break
    rows=[]
    for weight,cycle in sorted(selected,key=lambda item:(-item[0],item[1])):
        a,b,c,d=cycle;edges=[(a,b),(b,c),(c,d),(d,a)]
        for size in (1,3):
            for odd in itertools.combinations(range(4),size):
                terms=defaultdict(int)
                for index,(u,v) in enumerate(edges):
                    sign=1 if index in odd else -1
                    terms[u]+=sign;terms[v]+=sign
                    terms[products[tuple(sorted((u,v)))]]-=2*sign
                rows.append(({v:c for v,c in terms.items() if c},'<=',size-1))
    return rows,examined


def generate(model, upper, seconds, max_nodes, root_multipliers=None):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    rows=canonical_rows(model);variables=model.binary+list(model.products.values())
    indices={v:i for i,v in enumerate(variables)}
    groups=[[i for i,(_,sense,_) in enumerate(rows) if sense==target] for target in ('<=','=')]
    def matrix(group):
        triples=[(j,indices[v],c) for j,i in enumerate(group) for v,c in rows[i][0].items()]
        if not triples:return coo_matrix((len(group),len(variables))).tocsr()
        a,b,c=zip(*triples)
        return coo_matrix((c,(a,b)),shape=(len(group),len(variables))).tocsr()
    matrices=[matrix(g) for g in groups];rhs=[np.array([rows[i][2] for i in g],dtype=float) for g in groups]
    costs=np.array([model.objective.get(v,0) for v in variables],dtype=float)
    deadline=time.monotonic()+seconds;nodes=[];pending=[];solved=0;inherited_prunes=0
    def append(fixed,multipliers):
        index=len(nodes);lower=bound(model,rows,fixed,multipliers)
        nodes.append(dict(fixed=fixed,multipliers=list(multipliers),lower=lower))
        heapq.heappush(pending,(lower,index));return index
    append({},root_multipliers if root_multipliers is not None else ['0']*len(rows))
    while pending and time.monotonic()<deadline:
        _,index=heapq.heappop(pending);node=nodes[index]
        if node['lower']>=upper:inherited_prunes+=1;continue
        limits=[(node['fixed'].get(v,0),node['fixed'].get(v,1)) for v in variables]
        remaining=deadline-time.monotonic()
        if remaining<=0:break
        result=linprog(costs,A_ub=matrices[0] if groups[0] else None,b_ub=rhs[0] if groups[0] else None,
                       A_eq=matrices[1] if groups[1] else None,b_eq=rhs[1] if groups[1] else None,
                       bounds=limits,method='highs',options={'time_limit':min(10.,remaining),'threads':1})
        solved+=1
        # Infeasibility and floating status alone never close a proof obligation.
        if result.success:
            proposed=[Fraction(0)]*len(rows)
            for group,part in zip(groups,(result.ineqlin,result.eqlin)):
                for i,value in zip(group,part.marginals):proposed[i]=Fraction(float(-value)).limit_denominator(1000000)
            lower=bound(model,rows,node['fixed'],proposed)
            if lower>=node['lower']:node.update(multipliers=[str(x) for x in proposed],lower=lower)
        if node['lower']>=upper:continue
        free=[v for v in model.binary if v not in node['fixed']]
        if not free or len(nodes)+2>max_nodes or time.monotonic()>=deadline:continue
        variable=min(free,key=lambda v:(abs(result.x[indices[v]]-.5),v)) if result.success else free[0]
        node['branch']=variable
        # Parent proof is valid in both child domains. No previously established
        # bound is discarded when the tree representation changes.
        node['children']=[append({**node['fixed'],variable:bit},node['multipliers']) for bit in (0,1)]
    return dict(nodes=nodes,lp_calls=solved,inherited_prunes=inherited_prunes,deadline=time.monotonic()>=deadline,
                variables=len(variables),binary_variables=len(model.binary),rows=len(rows),strategy='weakest-bound-first; inherited dual witnesses')


def run(action,fixture,layout,certificate,seconds=120,max_nodes=101,cycles=500):
    evidence=Evidence(read_json(fixture));candidate=read_json(layout);upper=evidence.score(candidate)
    model=IntegerModel(evidence)
    if model.evaluate(model.assignment(candidate))!=upper:raise ValueError('Geometry/model disagreement')
    hashes={name:digest(path) for name,path in [('fixture_sha256',fixture),('layout_sha256',layout),
      ('model_checker_sha256',Path(__file__).with_name('proof.py')),('rational_checker_sha256',Path(__file__).with_name('order_certificate.py')),
      ('structural_checker_sha256',Path(__file__).with_name('structural.py')),('strengthened_checker_sha256',Path(__file__))]}
    saved=read_json(certificate) if action=='verify' else None
    if saved:
        if any(saved.get(k)!=v for k,v in hashes.items()):raise ValueError('Input/checker fingerprint mismatch')
        cycles=saved['cycle_limit']
    extra,examined=cycle_cuts(model,cycles);model.rows.extend(extra)
    base=structural_bound(evidence)['lower']
    if action=='check':saved={**hashes,**generate(model,upper,seconds,max_nodes), 'cycle_limit':cycles,'cycle_cuts':len(extra),'examined_cycles':examined}
    if saved['cycle_cuts']!=len(extra) or saved['examined_cycles']!=examined:raise ValueError('Cycle selection metadata differs')
    summary=replay(model,upper,saved);summary['ordering_lower']=summary['lower'];summary['structural_lower']=base
    summary['lower']=max(summary['lower'],base)
    if summary['lower']>upper:raise ValueError('Bound exceeds candidate')
    summary['gap']=upper-summary['lower'];summary['status']='EXACT_ORDER_CERTIFICATE' if summary['gap']==0 else 'UNRESOLVED_ORDER_GAP'
    if action=='check':saved['summary']=summary;certificate.write_text(json.dumps(saved,separators=(',',':'))+'\n')
    elif saved['summary']!=summary:raise ValueError('Conclusion differs from exact replay')
    print(('VERIFIED_' if action=='verify' else '')+summary['status'],f"U={upper} L={summary['lower']} gap={summary['gap']} nodes={summary['nodes']} cuts={len(extra)}",flush=True)
    return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['check','verify']);p.add_argument('fixture',type=Path);p.add_argument('layout',type=Path);p.add_argument('certificate',type=Path);p.add_argument('--seconds',type=float,default=120);p.add_argument('--max-nodes',type=int,default=101);p.add_argument('--cycles',type=int,default=500);a=p.parse_args()
    if a.seconds<0 or a.max_nodes<1 or a.cycles<0:p.error('Invalid budget')
    run(a.action,a.fixture,a.layout,a.certificate,a.seconds,a.max_nodes,a.cycles)
