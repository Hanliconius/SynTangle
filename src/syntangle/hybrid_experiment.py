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


def improve_neighborhoods(model, state, seconds, adaptive=False, progress=None, *, deduplicate=False):
    deadline=time.perf_counter()+seconds;records=[]
    species=model.fixture.species_ids
    for width in ((3,5) if adaptive else (3,)):
        while time.perf_counter()<deadline:
            changed=False;all_optimal=True;seen_windows=set()
            for center in range(len(species)):
                if time.perf_counter()>=deadline:break
                left=max(0,min(center-width//2,len(species)-width))
                window=tuple(species[left:left+width])
                if deduplicate and window in seen_windows:continue
                seen_windows.add(window)
                active=set(window)
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
                                    status=info['status'],fixed_decisions=len(fixed),
                                    nodes=info.get('nodes',0),solver_seconds=info.get('solver_seconds',0),
                                    api_preparation_seconds=info.get('api_preparation_seconds',0)))
                if progress and after<before:progress(state)
            # Stop repeating unchanged sweeps; adaptive version then expands.
            # This is a heuristic stopping decision, never a global exclusion.
            if not changed:break
    return state,dict(subsolves=len(records),neighborhoods=records,unique_windows_per_sweep=deduplicate,
                      stagnation_policy='expand 3 to 5 species' if adaptive else 'stop unchanged sweep')


def run_hybrid(fixture, starting_state, method, seconds, seed=1, progress=None, *, scheduling="weighted", backend="highs",
               mirror_symmetry=False, neighborhood_min_decisions=128, neighborhood_fraction=.25,
               neighborhood_max_seconds=60, deduplicate_neighborhoods=False):
    if method not in ('milp_reclaim','milp_hint','hybrid_3','hybrid_adaptive'):
        raise ValueError(method)
    if scheduling not in ("weighted", "equal") or backend not in ("highs", "scipy"):
        raise ValueError("Unknown scheduling/backend")
    if backend == "scipy" and method != "milp_reclaim":
        raise ValueError("SciPy comparison does not support a feasible MIP start")
    if not 0<=neighborhood_fraction<=1 or neighborhood_max_seconds<0 or neighborhood_min_decisions<0:
        raise ValueError('Invalid neighborhood limits')
    if mirror_symmetry and backend!='highs':raise ValueError('Mirror reduction requires the direct global API')
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
        model_started=time.perf_counter()
        model=JointModel(fixture,nodes,refs)
        model.graph_preparation_seconds=time.perf_counter()-model_started
        models.append((i,model))
    preparation=time.perf_counter()-started
    # Easy components first: their unused time flows directly to larger ones.
    models.sort(key=lambda item:(item[1].n+len(item[1].quadratic),item[0]))
    weights=[max(1,m.n+len(m.quadratic)) if scheduling=="weighted" else 1 for _,m in models]
    for position,(i,model) in enumerate(models):
        record=dict(component=i,decision_variables=model.n,crossing_factors=model.graph.factor_count,
            hard_orientation_groups=[[r.label for r in g] for g in model.basis.free_flip_groups],
            graph_preparation_seconds=model.graph_preparation_seconds)
        remaining=deadline-time.perf_counter()
        if remaining<=0:
            records.append({**record,'status':'deadline unstarted','component_lower':0});continue
        allocation=remaining*weights[position]/sum(weights[position:])
        component_deadline=time.perf_counter()+allocation
        component_start=time.perf_counter()
        # Small models already solve quickly: no costly heuristic prelude.
        if method.startswith('hybrid') and model.n>neighborhood_min_decisions:
            lns_started=time.perf_counter()
            state,info=improve_neighborhoods(model,state,min(neighborhood_max_seconds,allocation*neighborhood_fraction),
                                            method=='hybrid_adaptive',progress,
                                            deduplicate=deduplicate_neighborhoods)
            record.update(neighborhood_search=info,neighborhood_seconds=time.perf_counter()-lns_started)
        global_budget=max(.001,component_deadline-time.perf_counter())
        global_started=time.perf_counter()
        if backend=='highs':
            state,bound,info=model.highs(state,global_budget,progress=progress,
                                        use_mip_start=method!='milp_reclaim',mirror_symmetry=mirror_symmetry)
        else:
            state,bound,info=model.milp(state,global_budget)
        record['global_seconds']=time.perf_counter()-global_started
        lower+=bound;record.update(info,component_seconds=time.perf_counter()-component_start,
                                  allocated_seconds=allocation,global_budget_seconds=global_budget)
        records.append(record)
        if progress:progress(state)
    validate_saved_layout(fixture,state,whole)
    upper=score_crossings(fixture,state).crossings
    if upper>score_crossings(fixture,starting_state).crossings or lower>upper:
        raise AssertionError('Hybrid reconstruction/bounds lost incumbent')
    from .tangledness import tangledness_metrics
    metrics=tangledness_metrics(score_crossings(fixture,starting_state).crossings,upper,lower)
    return dict(method=method,optimized_state=state.to_dict(),starting_state=starting_state.to_dict(),
        tangledness=metrics,
        upper_bound=upper,lower_bound=lower,optimality_gap=upper-lower,
        optimality_status='proven optimum' if upper==lower else 'bounded best known',
        seconds=time.perf_counter()-started,preparation_seconds=preparation,
        scheduling=scheduling,backend=backend,
        requested_mirror_symmetry=mirror_symmetry,
        neighborhood_policy=dict(min_decisions=neighborhood_min_decisions,fraction=neighborhood_fraction,
            max_seconds=neighborhood_max_seconds,deduplicate=deduplicate_neighborhoods),
        component_diagnostics=records,
        whole_chromosome_order_changes={sp:[r.chromosome_id for r in state.chromosome_order[sp]]
            for sp in fixture.species_ids if state.chromosome_order[sp]!=starting_state.chromosome_order[sp]},
        whole_chromosome_flips=[r.label for r in fixture.chromosome_refs
            if state.chromosome_orientation[r]!=starting_state.chromosome_orientation[r]],
        reduction_policy='shared immutable graph/GF2/model; conditional neighborhoods then one global exact solve',
        global_search_restarts=0)
