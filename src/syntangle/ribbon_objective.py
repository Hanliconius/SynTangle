"""Ribbon-aware whole-chromosome refinement. Legacy midpoint objective is unchanged."""
import math,time
import numpy as np
from itertools import combinations
from .layout import score_crossings
from .orientation_space import orientation_basis
from .saved_layout import validate_saved_layout
from .visualize import _unambiguous_links


def twist_terms(fixture, model, width_weighted=False):
    links=[link for a,b in zip(fixture.species_ids,fixture.species_ids[1:]) for link in _unambiguous_links(fixture,a,b)]
    widths=[min(b1.end-b1.start,b2.end-b2.start) for _,(_,b1),(_,b2) in links]
    mean=sum(widths)/len(widths) if widths else 1
    group={r:i for i,g in enumerate(model.basis.free_flip_groups) for r in g}
    constant=0.;linear=np.zeros(model.n);quadratic={}
    for link,width in zip(links,widths):
        _,(c1,b1),(c2,b2)=link;a,b=group[c1.ref],group[c2.ref]
        base=(b1.strand=='-') ^ (b2.strand=='-') ^ (model.basis.base_assignment[c1.ref]<0) ^ (model.basis.base_assignment[c2.ref]<0)
        w=width/mean if width_weighted else 1.
        if a==b:constant+=w*base;continue
        a,b=sorted((a,b));constant+=w*base
        sign=-1 if base else 1
        linear[a]+=w*sign;linear[b]+=w*sign
        quadratic[a,b]=quadratic.get((a,b),0)-2*w*sign
    return constant,linear,quadratic


def ribbon_metrics(fixture,state):
    rows=[]
    for a,b in zip(fixture.species_ids,fixture.species_ids[1:]):
        for _,(c1,b1),(c2,b2) in _unambiguous_links(fixture,a,b):
            twist=(b1.strand=='-') ^ (b2.strand=='-') ^ (state.chromosome_orientation[c1.ref]<0) ^ (state.chromosome_orientation[c2.ref]<0)
            rows.append((twist,min(b1.end-b1.start,b2.end-b2.start)))
    mean=sum(w for _,w in rows)/len(rows) if rows else 1
    return dict(crossings=score_crossings(fixture,state).crossings,twists=sum(t for t,_ in rows),width_weighted_twists=sum(w/mean for t,w in rows if t))


def optimize_ribbons(fixture,start,mode='lex',weight=1.,width_weighted=False,seconds=600):
    if mode not in ('lex','weighted') or not math.isfinite(weight) or weight<0 or not math.isfinite(seconds) or seconds<=0:raise ValueError('Invalid ribbon refinement options')
    from .global_experiment import JointModel
    from .heuristic import _component_map
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix
    begin=time.monotonic();validate_saved_layout(fixture,start,orientation_basis(fixture,fixture.chromosome_refs))
    components,_=_component_map(fixture);nodes=tuple(n for comp in components for n in comp)
    model=JointModel(fixture,nodes,fixture.chromosome_refs)
    tc,tl,tq=twist_terms(fixture,model,width_weighted)
    pairs=sorted(set(model.quadratic)|set(tq));n=model.n
    crossings=np.r_[model.linear,[model.quadratic.get(p,0) for p in pairs]]
    twists=np.r_[tl,[tq.get(p,0) for p in pairs]]
    objective=twists if mode=='lex' else crossings+weight*twists
    objective_constant=tc if mode=='lex' else model.constant+weight*tc
    rr=[];cc=[];vv=[];lo=[];hi=[]
    def add(entries,lower,upper):
        row=len(lo);lo.append(lower);hi.append(upper)
        for col,val in entries:rr.append(row);cc.append(col);vv.append(val)
    for a,b,c in model.triangles:add([(a,1),(b,1),(c,-1)],0,1)
    for z,(a,b) in enumerate(pairs,n):
        add([(z,1),(a,-1)],-np.inf,0);add([(z,1),(b,-1)],-np.inf,0);add([(z,1),(a,-1),(b,-1)],-1,np.inf)
    cap=score_crossings(fixture,start).crossings
    if mode=='lex':add([(i,v) for i,v in enumerate(crossings) if v],-np.inf,cap-model.constant)
    matrix=coo_matrix((vv,(rr,cc)),shape=(len(lo),len(objective))).tocsr()
    bits=model.encode(start);initial_vec=np.r_[bits,[bits[a]*bits[b] for a,b in pairs]]
    result=milp(objective,integrality=np.r_[np.ones(n),np.zeros(len(pairs))],bounds=Bounds(np.zeros(len(objective)),np.ones(len(objective))),constraints=LinearConstraint(matrix,lo,hi),options=dict(time_limit=max(.001,seconds-(time.monotonic()-begin)),mip_rel_gap=0,threads=1))
    final=start
    if result.x is not None:
        candidate=model.decode(np.rint(result.x[:n]),start)
        cbits=model.encode(candidate);vec=np.r_[cbits,[cbits[a]*cbits[b] for a,b in pairs]]
        if objective@vec<=objective@initial_vec+1e-7:final=candidate
    before=ribbon_metrics(fixture,start);after=ribbon_metrics(fixture,final)
    if mode=='lex' and after['crossings']>cap:raise AssertionError('Crossing cap violated')
    validate_saved_layout(fixture,final,orientation_basis(fixture,fixture.chromosome_refs))
    return dict(mode=mode,weight=weight,width_weighted=width_weighted,starting_metrics=before,final_metrics=after,optimized_state=final.to_dict(),crossing_cap=cap if mode=='lex' else None,solver_status=int(result.status),solver_message=result.message,objective_constant=objective_constant,objective_dual_bound=float(result.mip_dual_bound+objective_constant) if getattr(result,'mip_dual_bound',None) is not None else None,optimality_status='numerical solver result; no independent certificate',seconds=time.monotonic()-begin)


def refine_ribbons(args):
    """Public file interface; records both objectives and immutable input hashes."""
    import csv
    import hashlib
    import json
    from pathlib import Path
    from .fixtures import load_fixture
    from .genespace import import_layout
    from .saved_layout import decode_saved_layout

    fixture = load_fixture(args.fixture)
    start = import_layout(fixture, args.layout)
    result = optimize_ribbons(
        fixture, start, mode='lex' if args.objective == 'crossing-cap' else 'weighted',
        weight=args.weight, width_weighted=args.twist_metric == 'width', seconds=args.seconds)
    final = decode_saved_layout(fixture, result['optimized_state'])
    result.update(objective=args.objective, twist_metric=args.twist_metric,
                  fixture_sha256=hashlib.sha256(Path(args.fixture).read_bytes()).hexdigest(),
                  layout_input_sha256=hashlib.sha256(Path(args.layout).read_bytes()).hexdigest(),
                  allowed_moves='whole-chromosome permutations and flips; fixed internal coordinates and links')
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'layout.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['species', 'chromosome', 'rank', 'orientation'])
        for species, order in final.chromosome_order.items():
            for rank, ref in enumerate(order, 1):
                writer.writerow([species, ref.chromosome_id, rank, final.chromosome_orientation[ref]])
    (out / 'audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(before=result['starting_metrics'], after=result['final_metrics'],
                          objective=args.objective, output=str(out))))
    return 0
