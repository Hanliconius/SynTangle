"""Root separation of valid four-cycle inequalities with exact rational replay.

The LP solution chooses useful cuts, never which chromosome layouts are legal.
Every selected row is regenerated from the original integer model by the verifier.
"""
import argparse,json,time
from fractions import Fraction
from pathlib import Path
from proof import Evidence,IntegerModel,read_json,digest
from order_certificate import canonical_rows,bound,replay
from order_strengthened import cycle_cuts


def solve(model,seconds):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    rows=canonical_rows(model);variables=model.binary+list(model.products.values());indices={v:i for i,v in enumerate(variables)}
    groups=[[i for i,r in enumerate(rows) if r[1]==sense] for sense in ('<=','=')]
    def matrix(group):
        triples=[(j,indices[v],c) for j,i in enumerate(group) for v,c in rows[i][0].items()]
        if not triples:return coo_matrix((len(group),len(variables))).tocsr()
        a,b,c=zip(*triples);return coo_matrix((c,(a,b)),shape=(len(group),len(variables))).tocsr()
    result=linprog([model.objective.get(v,0) for v in variables],
      A_ub=matrix(groups[0]) if groups[0] else None,b_ub=[rows[i][2] for i in groups[0]] if groups[0] else None,
      A_eq=matrix(groups[1]) if groups[1] else None,b_eq=[rows[i][2] for i in groups[1]] if groups[1] else None,
      bounds=(0,1),method='highs',options={'time_limit':max(.001,seconds),'threads':1})
    weights=[Fraction(0)]*len(rows)
    if result.success:
        for group,part in zip(groups,(result.ineqlin,result.eqlin)):
            for i,value in zip(group,part.marginals):weights[i]=Fraction(float(-value)).limit_denominator(1000000)
    return rows,weights,({v:float(result.x[i]) for v,i in indices.items()} if result.success else None)


def violations(pool,values,selected,batch):
    used=set(selected)
    ranked=[(sum(c*values[v] for v,c in terms.items())-rhs,i) for i,(terms,_,rhs) in enumerate(pool) if i not in used]
    return [i for value,i in sorted(ranked,key=lambda x:(-x[0],x[1]))[:batch] if value>1e-7]


def run(action,fixture,layout,certificate,seconds=60,batch=32,rounds=20,cycles=500):
    evidence=Evidence(read_json(fixture));upper=evidence.score(read_json(layout));model=IntegerModel(evidence)
    hashes={str(n):digest(p) for n,p in [('fixture_sha256',fixture),('layout_sha256',layout),('selective_checker_sha256',Path(__file__)),
      ('model_checker_sha256',Path(__file__).with_name('proof.py')),('rational_checker_sha256',Path(__file__).with_name('order_certificate.py')),('cycle_checker_sha256',Path(__file__).with_name('order_strengthened.py'))]}
    saved=read_json(certificate) if action=='verify' else None
    if saved:
        if any(saved.get(k)!=v for k,v in hashes.items()):raise ValueError('Input/checker changed')
        cycles=saved['cycle_limit']
    pool,examined=cycle_cuts(model,cycles)
    if action=='check':
        deadline=time.monotonic()+seconds;selected=[];history=[];rows=canonical_rows(model);weights=['0']*len(rows);lower=bound(model,rows,{},weights)
        for step in range(rounds):
            if time.monotonic()>=deadline:break
            rows,proposed,values=solve(model,min(10,deadline-time.monotonic()))
            candidate=bound(model,rows,{},proposed)
            if candidate>=lower:lower=candidate;weights=[str(x) for x in proposed]
            history.append(dict(round=step,lower=lower,cuts=len(selected)))
            if lower>=upper or values is None or step+1==rounds:break
            added=violations(pool,values,selected,batch)
            if not added:break
            selected.extend(added);model.rows.extend(pool[i] for i in added)
            weights.extend(['0']*len(added)) # previous exact witness survives every added row
        saved={**hashes,'selected':selected,'cycle_limit':cycles,'examined_cycles':examined,'pool_cuts':len(pool),'history':history,
          'nodes':[dict(fixed={},multipliers=weights,lower=lower)],'lp_calls':len(history),'cycle_cuts':len(selected),'rows':len(canonical_rows(model))}
    else:
        selected=saved['selected']
        if len(set(selected))!=len(selected) or any(type(i) is not int or not 0<=i<len(pool) for i in selected):raise ValueError('Invalid selected cut indices')
        if saved['pool_cuts']!=len(pool) or saved['examined_cycles']!=examined or saved['cycle_cuts']!=len(selected):raise ValueError('Cut pool mismatch')
        model.rows.extend(pool[i] for i in selected)
    summary=replay(model,upper,saved)
    if action=='check':saved['summary']=summary;certificate.write_text(json.dumps(saved,separators=(',',':'))+'\n')
    elif summary!=saved['summary']:raise ValueError('Conclusion differs')
    print(('VERIFIED_' if action=='verify' else '')+summary['status'],f"U={upper} L={summary['lower']} gap={summary['gap']} nodes=1 cuts={len(selected)} LP={saved['lp_calls']}",flush=True)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['check','verify']);p.add_argument('fixture',type=Path);p.add_argument('layout',type=Path);p.add_argument('certificate',type=Path)
    p.add_argument('--seconds',type=float,default=60);p.add_argument('--batch',type=int,default=32);p.add_argument('--rounds',type=int,default=20);p.add_argument('--cycles',type=int,default=500);a=p.parse_args()
    if a.seconds<0 or min(a.batch,a.rounds)<1 or a.cycles<0:p.error('Invalid budget')
    run(a.action,a.fixture,a.layout,a.certificate,a.seconds,a.batch,a.rounds,a.cycles)
