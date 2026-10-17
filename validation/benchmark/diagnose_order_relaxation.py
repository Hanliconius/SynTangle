"""Conditional gap probes, not bounds for an unconstrained whole-fixture solve."""
import argparse
import csv
from itertools import combinations, permutations
import json
from pathlib import Path
import time


def pair_minimum(costs):
    return sum(min(costs.cost(a,b), costs.cost(b,a)) for a,b in combinations(costs.refs,2))


def triangle_packing(costs):
    """Admissible transitivity correction using edge-disjoint triangles."""
    refs = costs.refs
    index = {ref: i for i,ref in enumerate(refs)}
    candidates = []
    for triple in combinations(refs,3):
        baseline = sum(min(costs.cost(a,b), costs.cost(b,a)) for a,b in combinations(triple,2))
        legal = min(sum(costs.cost(a,b) for i,a in enumerate(order) for b in order[i+1:])
                    for order in permutations(triple))
        gain = legal - baseline
        if gain:
            edges = frozenset(tuple(sorted((index[a],index[b]))) for a,b in combinations(triple,2))
            candidates.append((gain, tuple(index[r] for r in triple), edges))
    used = set(); correction = 0; packed = 0
    for gain, _, edges in sorted(candidates, key=lambda row: (-row[0], row[1])):
        if not edges & used:
            correction += gain
            packed += 1
            used.update(edges)
    return correction, packed, len(candidates)


def probe_row(fixture, state, nodes, species):
    from syntangle.incidence import chromosome_node_id
    from syntangle.order_dp import build_pairwise_order_costs, solve_order_subset_dp
    from syntangle.pair_cost import pair_component_crossings
    index = fixture.species_ids.index(species)
    neighbors = tuple(fixture.species_ids[i] for i in (index-1,index+1)
                      if 0 <= i < len(fixture.species_ids))
    current = tuple(r for r in state.chromosome_order[species]
                    if chromosome_node_id(r) in nodes)
    def actual(neighbor):
        other = tuple(r for r in state.chromosome_order[neighbor]
                      if chromosome_node_id(r) in nodes)
        if fixture.species_ids.index(neighbor) < index:
            return pair_component_crossings(fixture, neighbor, species, other, current,
                                            state.chromosome_orientation, nodes)
        return pair_component_crossings(fixture, species, neighbor, current, other,
                                        state.chromosome_orientation, nodes)
    current_cost = sum(actual(n) for n in neighbors)
    singles = [build_pairwise_order_costs(fixture,state,species,nodes,neighbors=(n,)) for n in neighbors]
    combined = build_pairwise_order_costs(fixture,state,species,nodes,neighbors=neighbors)
    current_variable = sum(combined.cost(a,b) for i,a in enumerate(current) for b in current[i+1:])
    constant = current_cost - current_variable
    if constant < 0:
        raise AssertionError('Conditional pair decomposition disagrees with canonical edge scoring')
    separate = constant + sum(pair_minimum(costs) for costs in singles)
    joint = constant + pair_minimum(combined)
    correction, packed, conflicting = triangle_packing(combined)
    triangle = joint + correction
    exact = None
    if len(current) <= 12:
        exact = constant + solve_order_subset_dp(species, combined).variable_crossing_cost
    if not 0 <= separate <= joint <= triangle <= current_cost:
        raise AssertionError('Invalid conditional bound hierarchy')
    if exact is not None and not triangle <= exact <= current_cost:
        raise AssertionError('Triangle bound exceeds exact conditional optimum')
    return dict(species=species, chromosomes=len(current), current_adjacent_crossings=current_cost,
                separate_neighbor_pair_bound=separate, joint_neighbor_pair_bound=joint,
                neighbor_consistency_lift=joint-separate,
                triangle_bound=triangle, transitivity_lift=correction,
                edge_disjoint_triangles=packed, positive_triangle_count=conflicting,
                exact_conditional_optimum=exact,
                note='Orientations and neighboring orders fixed; not a global lower bound')


def main():
    from syntangle import load_validation_bundle
    from syntangle.saved_layout import decode_saved_layout
    from syntangle.heuristic import _component_map
    from syntangle.bounds import build_relaxed_crossing_bound
    from syntangle.orientation_space import orientation_basis
    from benchmark_bound_pair import choose_saved
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',required=True);parser.add_argument('--history',required=True)
    parser.add_argument('--output',required=True);parser.add_argument('--index',type=int,required=True)
    args=parser.parse_args();root=Path(args.root)
    with (root/'benchmark_manifest.tsv').open() as handle:
        entry=list(csv.DictReader(handle,delimiter='\t'))[args.index]
    fixture=load_validation_bundle(root/entry['case_dir'])
    score,source,data=choose_saved(fixture,entry['case_id'],Path(args.history))
    state=decode_saved_layout(fixture,data)
    components,component_of=_component_map(fixture)
    report=dict(case_id=entry['case_id'],saved_crossings=score,source=source,components=[])
    prior=json.loads(Path(source).read_text())
    candidates=[prior]+[row.get('result',{}) for row in prior.get('methods',[])]
    records=[record for record in candidates
             if record.get('optimized_score',{}).get('crossings')==score]
    if records:
        diagnostics=records[0].get('component_diagnostics',[])
        counters=dict(table_entries=sum(d.get('table_entries_evaluated',0) for d in diagnostics),
                      branch_nodes=sum(d.get('nodes_evaluated',0) for d in diagnostics
                                       if d.get('method')=='implicit-branch-and-bound'),
                      continuation_nodes=sum(d.get('continuation_nodes',0) for d in diagnostics))
        report['prior_result_counters']=counters
        print('PRIOR_RESULT_WORK '+ ' '.join(f'{key}={value}' for key,value in counters.items()),flush=True)

    print(f"{entry['case_id']} SAVED_C={score}; all row probes are conditional",flush=True)
    for i,nodes in enumerate(components):
        refs=tuple(sorted(r for r in fixture.chromosome_refs if component_of[r]==i))
        basis=orientation_basis(fixture,refs)
        bound=build_relaxed_crossing_bound(fixture,nodes,basis)
        bits=tuple(int(state.chromosome_orientation[group[0]] != basis.base_assignment[group[0]])
                   for group in basis.free_flip_groups)
        component=dict(component_id=i,root_relaxed_bound=bound.lower_bound((None,)*len(bits)),
                       saved_orientation_relaxed_bound=bound.lower_bound(bits),rows=[])
        if len(refs) > len(fixture.species_ids):
            for species in fixture.species_ids:
                if not any(r.species_id == species for r in refs):
                    continue
                started=time.perf_counter();row=probe_row(fixture,state,nodes,species)
                row['seconds']=time.perf_counter()-started;component['rows'].append(row)
                print(f"component={i} {species} n={row['chromosomes']} current={row['current_adjacent_crossings']} "
                      f"separate={row['separate_neighbor_pair_bound']} joint={row['joint_neighbor_pair_bound']} "
                      f"triangle={row['triangle_bound']} exact={row['exact_conditional_optimum']} "
                      f"neighbor_lift={row['neighbor_consistency_lift']} triangle_lift={row['transitivity_lift']} "
                      f"seconds={row['seconds']:.3f}",flush=True)
        report['components'].append(component)
        print(f"component={i} root_orientation_bound={component['root_relaxed_bound']} "
              f"saved_orientation_bound={component['saved_orientation_relaxed_bound']}",flush=True)
        output=Path(args.output)/entry['case_id'];output.mkdir(parents=True,exist_ok=True)
        (output/'diagnostic.json').write_text(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
