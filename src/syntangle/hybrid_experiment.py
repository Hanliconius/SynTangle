"""Graph preprocessing -> optional conditional neighborhoods -> one global solve.

Global exact search starts once per component; no search frontier is restarted.
Neighborhood exclusions are scoped to fixed outside decisions, never carried as
unconditional restrictions. Immutable crossing/MILP matrices are shared.
"""
import time
from .global_experiment import JointModel
from .heuristic import _component_map
from .orientation_space import orientation_basis
from .saved_layout import validate_saved_layout
from .layout import canonicalize_component_order, score_crossings


def improve_neighborhoods(model, state, seconds, adaptive=False, progress=None):
    deadline=time.perf_counter()+seconds;records=[]
    species=model.fixture.species_ids
    for width in ((3,5) if adaptive else (3,)):
        while time.perf_counter()<deadline:
            changed=False;all_optimal=True
            for center in range(len(species)):
                if time.perf_counter()>=deadline:break
                left=max(0,min(center-width//2,len(species)-width))
                active=set(species[left:left+width])
                bits=model.encode(state)
                free={j for (a,b),j in model.orders.items() if a.species_id in active}
                free.update(j for j,group in enumerate(model.basis.free_flip_groups)
                            if any(r.species_id in active for r in group))
                fixed={j:value for j,value in enumerate(bits) if j not in free}
                before=model.objective(bits)
                state,_,info=model.milp(state,min(5,deadline-time.perf_counter()),fixed)
                after=model.objective(model.encode(state))
                changed |= after<before;all_optimal &= info['status']==0
                records.append(dict(width=width,center=center,before=int(before),after=int(after),
                                    status=info['status'],fixed_decisions=len(fixed)))
                if progress and after<before:progress(state)
            # Stop repeating unchanged sweeps; adaptive version then expands.
            # This is a heuristic stopping decision, never a global exclusion.
            if not changed:break
    return state,dict(subsolves=len(records),neighborhoods=records,
                      stagnation_policy='expand 3 to 5 species' if adaptive else 'stop unchanged sweep')


def run_hybrid(fixture, starting_state, method, seconds, seed=1, progress=None):
    if method not in ('milp_reclaim','milp_hint','hybrid_3','hybrid_adaptive'):
        raise ValueError(method)
    started=time.perf_counter();deadline=started+seconds
    whole=orientation_basis(fixture,fixture.chromosome_refs)
    validate_saved_layout(fixture,starting_state,whole)
    state=canonicalize_component_order(fixture,starting_state)
    components,component_of=_component_map(fixture)
    models=[];records=[];lower=0
    for i,nodes in enumerate(components):
        if time.perf_counter()>=deadline:
            records.extend(dict(component=j,status='preparation deadline unstarted',component_lower=0)
                           for j in range(i,len(components)))
            break
        refs=tuple(r for r in fixture.chromosome_refs if component_of[r]==i)
        if score_crossings(fixture,state,restrict_component_nodes=nodes).crossings==0:
            records.append(dict(component=i,status='proven zero; retained and removed from search',
                component_lower=0,component_upper=0,retained_chromosomes=[r.label for r in refs],
                reason='legal zero-crossing incumbent and nonnegative crossing objective'))
            continue
        models.append((i,JointModel(fixture,nodes,refs)))
    preparation=time.perf_counter()-started
    # Easy components first: their unused time flows directly to larger ones.
    models.sort(key=lambda item:(item[1].n+len(item[1].quadratic),item[0]))
    weights=[max(1,m.n+len(m.quadratic)) for _,m in models]
    for position,(i,model) in enumerate(models):
        record=dict(component=i,decision_variables=model.n,crossing_factors=model.graph.factor_count,
            hard_orientation_groups=[[r.label for r in g] for g in model.basis.free_flip_groups])
        remaining=deadline-time.perf_counter()
        if remaining<=0:
            records.append({**record,'status':'deadline unstarted','component_lower':0});continue
        allocation=remaining*weights[position]/sum(weights[position:])
        component_deadline=time.perf_counter()+allocation
        component_start=time.perf_counter()
        # Small models already solve quickly: no costly heuristic prelude.
        if method.startswith('hybrid') and model.n>128:
            lns_started=time.perf_counter()
            state,info=improve_neighborhoods(model,state,min(60,allocation*.25),
                                            method=='hybrid_adaptive',progress)
            record.update(neighborhood_search=info,neighborhood_seconds=time.perf_counter()-lns_started)
        global_budget=max(.001,component_deadline-time.perf_counter())
        state,bound,info=model.highs(state,global_budget,progress=progress,
                                    use_mip_start=method!='milp_reclaim')
        lower+=bound;record.update(info,component_seconds=time.perf_counter()-component_start,
                                  allocated_seconds=allocation,global_budget_seconds=global_budget)
        records.append(record)
        if progress:progress(state)
    validate_saved_layout(fixture,state,whole)
    upper=score_crossings(fixture,state).crossings
    if upper>score_crossings(fixture,starting_state).crossings or lower>upper:
        raise AssertionError('Hybrid reconstruction/bounds lost incumbent')
    return dict(method=method,optimized_state=state.to_dict(),starting_state=starting_state.to_dict(),
        upper_bound=upper,lower_bound=lower,optimality_gap=upper-lower,
        optimality_status='proven optimum' if upper==lower else 'bounded best known',
        seconds=time.perf_counter()-started,preparation_seconds=preparation,
        component_diagnostics=records,
        whole_chromosome_order_changes={sp:[r.chromosome_id for r in state.chromosome_order[sp]]
            for sp in fixture.species_ids if state.chromosome_order[sp]!=starting_state.chromosome_order[sp]},
        whole_chromosome_flips=[r.label for r in fixture.chromosome_refs
            if state.chromosome_orientation[r]!=starting_state.chromosome_orientation[r]],
        reduction_policy='shared immutable graph/GF2/model; conditional neighborhoods then one global exact solve',
        global_search_restarts=0)
