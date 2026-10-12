"""Continue a replayed selective certificate without discarding its tree or cuts."""
import argparse,copy,heapq,json,time
from fractions import Fraction
from pathlib import Path
from proof import Evidence,IntegerModel,read_json,digest
from order_certificate import bound,canonical_rows,replay
from order_strengthened import cycle_cuts
from order_selective import run as selective


def continue_tree(model,upper,prior,seconds,max_nodes):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    nodes=copy.deepcopy(prior['nodes']);rows=canonical_rows(model)
    variables=model.binary+list(model.products.values());indices={v:i for i,v in enumerate(variables)}
    groups=[[i for i,r in enumerate(rows) if r[1]==sense] for sense in ('<=','=')]
    def matrix(group):
        triples=[(j,indices[v],c) for j,i in enumerate(group) for v,c in rows[i][0].items()]
        if not triples:return coo_matrix((len(group),len(variables))).tocsr()
        a,b,c=zip(*triples);return coo_matrix((c,(a,b)),shape=(len(group),len(variables))).tocsr()
    matrices=[matrix(g) for g in groups];rhs=[[rows[i][2] for i in g] for g in groups]
    costs=[model.objective.get(v,0) for v in variables]
    # All open proof obligations survive. Certified closed leaves stay closed.
    pending=[(n['lower'],i) for i,n in enumerate(nodes) if 'children' not in n and n['lower']<upper]
    heapq.heapify(pending);deadline=time.monotonic()+seconds;calls=0
    while pending and time.monotonic()<deadline:
        _,index=heapq.heappop(pending);node=nodes[index]
        free=[v for v in model.binary if v not in node['fixed']]
        if not free or len(nodes)+2>max_nodes:continue
        remaining=deadline-time.monotonic()
        if remaining<=0:break
        limits=[(node['fixed'].get(v,0),node['fixed'].get(v,1)) for v in variables]
        result=linprog(costs,A_ub=matrices[0] if groups[0] else None,b_ub=rhs[0] if groups[0] else None,
            A_eq=matrices[1] if groups[1] else None,b_eq=rhs[1] if groups[1] else None,
            bounds=limits,method='highs',options={'time_limit':min(10.,remaining),'threads':1})
        calls+=1
        if result.success:
            proposed=[Fraction(0)]*len(rows)
            for group,part in zip(groups,(result.ineqlin,result.eqlin)):
                for i,value in zip(group,part.marginals):proposed[i]=Fraction(float(-value)).limit_denominator(1000000)
            lower=bound(model,rows,node['fixed'],proposed)
            if lower>=node['lower']:node.update(lower=lower,multipliers=[str(x) for x in proposed])
        if node['lower']>=upper or time.monotonic()>=deadline:continue
        variable=min(free,key=lambda v:(abs(result.x[indices[v]]-.5),v)) if result.success else free[0]
        node['branch']=variable;node['children']=[]
        for bit in (0,1):
            fixed={**node['fixed'],variable:bit};weights=list(node['multipliers']);lower=bound(model,rows,fixed,weights)
            child=len(nodes);nodes.append(dict(fixed=fixed,multipliers=weights,lower=lower));node['children'].append(child)
            if lower<upper:heapq.heappush(pending,(lower,child))
    return nodes,calls,time.monotonic()>=deadline


def run(fixture,layout,prior_path,out,seconds=900,max_nodes=10001):
    old=selective('verify',fixture,layout,prior_path)
    prior=read_json(prior_path);e=Evidence(read_json(fixture));model=IntegerModel(e)
    pool,_=cycle_cuts(model,prior['cycle_limit']);model.rows.extend(pool[i] for i in prior['selected'])
    nodes,calls,deadline=continue_tree(model,old['upper'],prior,seconds,max_nodes)
    saved=copy.deepcopy(prior);saved.update(nodes=nodes,lp_calls=prior['lp_calls']+calls,deadline=deadline,
        resumed_from_sha256=digest(prior_path),prior_lower=old['lower'],prior_nodes=len(prior['nodes']),
        additional_lp_calls=calls,resume_seconds=seconds,resume_checker_sha256=digest(Path(__file__)))
    saved['summary']=replay(model,old['upper'],saved)
    if saved['summary']['lower']<old['lower']:raise ValueError('Previously certified bound was lost')
    out.write_text(json.dumps(saved,separators=(',',':'))+'\n')
    final=selective('verify',fixture,layout,out)
    print('CONTINUATION',json.dumps(dict(**final,prior_lower=old['lower'],additional_lp_calls=calls)),flush=True)
    return final

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('fixture','layout','prior','output'):p.add_argument(n,type=Path)
    p.add_argument('--seconds',type=float,default=900);p.add_argument('--max-nodes',type=int,default=10001);a=p.parse_args()
    if a.seconds<0 or a.max_nodes<1:p.error('Invalid budgets')
    run(a.fixture,a.layout,a.prior,a.output,a.seconds,a.max_nodes)
