"""Monotone component dispatch and in-place reduced-factor continuation.

There is deliberately no whole-fixture exception/fallback restart. Materialization
is preflighted before table reductions, and elimination-budget exhaustion is
handled inside the reduced factor solver with its reconstruction stack intact.
"""
from .search_control import controlled, expired, _progress
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
from statistics import median

from .layout import (LayoutState, ExactLayoutResult, SearchSpaceTooLarge,
                     initial_layout_state, canonicalize_component_order, score_crossings)
from .orientation_space import OrientationBasis, orientation_basis
from .residual import build_residual_factorization
from .residual_solver import _solve_incidence_component
from .branch_bound import _solve_branch_component


def median_order_seed(fixture, normalized, components, component_of, basis):
    """GENESPACE-style median homology-position sweep; not native GENESPACE.

    Only whole chromosomes move. Hard orientation equations use the shared basis.
    Missing neighbors leave the chromosome in its current relative position.
    """
    state = LayoutState(dict(normalized.chromosome_order), dict(basis.base_assignment))
    chromosomes = {chrom.ref: chrom for chrom in fixture.chromosomes}
    species_ids = fixture.species_ids
    best = state
    best_score = score_crossings(fixture, state).crossings
    for species_sequence in (species_ids, tuple(reversed(species_ids))):
        for species in species_sequence:
            positions = {}
            index = species_ids.index(species)
            neighbors = species_ids[max(0, index - 1):index] + species_ids[index + 1:index + 2]
            for neighbor in neighbors:
                for rank, ref in enumerate(state.chromosome_order[neighbor]):
                    chrom = chromosomes[ref]
                    for block in chrom.blocks:
                        fraction = (block.start + block.end) / (2 * chrom.length)
                        if state.chromosome_orientation[ref] == -1:
                            fraction = 1 - fraction
                        positions.setdefault(block.homology_id, []).append(rank + fraction)
            order = dict(state.chromosome_order)
            rebuilt = []
            for component_id in range(len(components)):
                refs = [r for r in order[species] if component_of[r] == component_id]
                ranks = {r: i for i, r in enumerate(order[species])}
                def key(ref):
                    values = [x for b in chromosomes[ref].blocks
                              for x in positions.get(b.homology_id, ())]
                    return (median(values) if values else ranks[ref], ranks[ref], ref.label)
                rebuilt.extend(sorted(refs, key=key))
            order[species] = tuple(rebuilt)
            state = LayoutState(order, dict(state.chromosome_orientation))
            score = score_crossings(fixture, state).crossings
            if score < best_score:
                best, best_score = state, score
    return best


def _solve_component(args):
    (fixture, normalized, incumbent, component_id, nodes, component_of, basis,
     variables, factors, permutation_cap, transition_cap, node_cap) = args
    try:
        result = _solve_incidence_component(
            fixture, component_id, nodes, component_of, variables, factors,
            permutation_cap, transition_cap, transition_cap, 8,
            prepared_basis=basis, incumbent_state=incumbent,
            continuation_node_cap=node_cap)
        return dict(component_id=component_id, orders=result.orders,
                    orientation=result.orientation, upper=result.optimum,
                    lower=result.lower_bound,
                    nodes=result.diagnostics.table_entries_evaluated + result.continuation_nodes,
                    method='reduced-factor-continuation' if result.handoffs else 'residual-elimination',
                    diagnostic={**result.diagnostics.to_dict(),
                                'handoffs': list(result.handoffs),
                                'continuation_nodes': result.continuation_nodes,
                                'continuation_branches_pruned': result.continuation_pruned})
    except SearchSpaceTooLarge as exc:
        # With continuation enabled the only escaping cap failures occur before
        # any factor-table scope reduction or elimination. Never catch a later
        # failure and restart the original search.
        reason = str(exc)
        if not (reason.startswith('Residual order variable') or
                reason.startswith('Factor materialization exceeds cap')):
            raise AssertionError('Unsafe fallback after reductions') from exc
        result = _solve_branch_component(fixture, normalized, incumbent, component_id,
                                         nodes, component_of, node_cap, basis)
        return dict(component_id=component_id, orders=result.orders,
                    orientation=result.orientation, upper=result.diagnostic.upper_bound,
                    lower=result.lower_bound, nodes=result.nodes_evaluated,
                    method='implicit-branch-and-bound',
                    diagnostic={**result.diagnostic.to_dict(), 'preflight_reason': reason,
                                'factor_reductions_before_dispatch': 0})


@controlled
def optimize_pipeline(fixture, *, orientation_cap_per_component,
                      permutation_cap_per_species, transition_cap_per_component,
                      branch_node_cap_per_component, local_restarts,
                      local_max_improving_steps, seed, component_workers, progress_callback, time_limit_seconds=None):
    from .heuristic import AutoLayoutResult, _component_map, optimize_local_search
    for value in (orientation_cap_per_component, permutation_cap_per_species,
                  transition_cap_per_component, branch_node_cap_per_component,
                  component_workers, local_restarts):
        if value < 1:
            raise ValueError('Solver caps, workers and restarts must be at least 1')
    if time_limit_seconds is not None and component_workers != 1:
        raise ValueError('Deadline-controlled solves currently require component_workers=1')
    components, component_of = _component_map(fixture)
    bases = {i: orientation_basis(fixture, tuple(sorted(
        r for r in fixture.chromosome_refs if component_of[r] == i)))
        for i in range(len(components))}
    whole_basis = OrientationBasis(
        {r: sign for b in bases.values() for r, sign in b.base_assignment.items()},
        tuple(g for b in bases.values() for g in b.free_flip_groups))
    initial = initial_layout_state(fixture)
    normalized = canonicalize_component_order(fixture, initial)
    seed_state = median_order_seed(fixture, normalized, components, component_of, whole_basis)
    incumbent = optimize_local_search(
        fixture, restarts=local_restarts,
        # Keep the first pass bounded in iterations and exact subset-DP size.
        max_improving_steps=min(local_max_improving_steps, 20), seed=seed,
        order_dp_max_chromosomes=12, starting_state=seed_state,
        prepared_components=(components, component_of), prepared_basis=whole_basis,
        progress_callback=progress_callback).layout.optimized_state
    factorization = build_residual_factorization(
        fixture, prepared_components=(components, component_of), prepared_bases=bases)
    arguments = [(fixture, normalized, incumbent, i, nodes, component_of, bases[i],
                  tuple(v for v in factorization.variables if v.component_id == i),
                  tuple(f for f in factorization.factors if f.component_id == i),
                  permutation_cap_per_species, transition_cap_per_component,
                  branch_node_cap_per_component) for i, nodes in enumerate(components)]
    workers = min(component_workers, max(1, len(components)))
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn')) as pool:
            results = list(pool.map(_solve_component, arguments))
    else:
        results = []
        current = incumbent
        current_score = score_crossings(fixture, current).crossings
        for args in arguments:
            component_id = args[3]
            def checkpoint(candidate):
                nonlocal current, current_score
                order = {sp: tuple(r for r in candidate.chromosome_order[sp]
                                  if component_of[r] == component_id)
                         for sp in fixture.species_ids}
                rebuilt = {sp: tuple(r for i in range(len(components))
                                    for r in (order[sp] if i == component_id else
                                              tuple(x for x in current.chromosome_order[sp]
                                                    if component_of[x] == i)))
                           for sp in fixture.species_ids}
                signs = dict(current.chromosome_orientation)
                signs.update({r: candidate.chromosome_orientation[r]
                              for r in bases[component_id].base_assignment})
                combined = LayoutState(rebuilt, signs)
                value = score_crossings(fixture, combined).crossings
                if value < current_score:
                    current, current_score = combined, value
                    if progress_callback:
                        progress_callback(combined, value)
            if expired():
                # Unstarted components retain their legal incumbent and the
                # universal nonnegative bound. No search is restarted.
                result = dict(component_id=component_id,
                    orders={sp: tuple(r for r in current.chromosome_order[sp]
                                      if component_of[r] == component_id)
                            for sp in fixture.species_ids},
                    orientation={r: current.chromosome_orientation[r]
                                 for r in bases[component_id].base_assignment},
                    upper=score_crossings(fixture, current,
                        restrict_component_nodes=args[4]).crossings,
                    lower=0, nodes=0, method='deadline-unstarted',
                    diagnostic={'deadline_reached': True, 'search_started': False})
            else:
                token = _progress.set(checkpoint)
                try:
                    result = _solve_component(args)
                finally:
                    _progress.reset(token)
                candidate = LayoutState(
                    {sp: tuple(r for i in range(len(components)) for r in
                         (result['orders'].get(sp, ()) if i == component_id else
                          tuple(x for x in current.chromosome_order[sp] if component_of[x] == i)))
                     for sp in fixture.species_ids},
                    {**current.chromosome_orientation, **result['orientation']})
                checkpoint(candidate)
            results.append(result)
    orders = {species: tuple(ref for result in results
                             for ref in result['orders'].get(species, ()))
              for species in fixture.species_ids}
    optimized = LayoutState(orders, {ref: sign for result in results
                                    for ref, sign in result['orientation'].items()})
    score = score_crossings(fixture, optimized)
    upper = sum(result['upper'] for result in results)
    lower = sum(result['lower'] for result in results)
    if score.crossings != upper or not 0 <= lower <= upper:
        raise AssertionError('Reduced component reconstruction/bounds disagree with canonical scoring')
    if upper > score_crossings(fixture, incumbent).crossings:
        raise AssertionError('Continuation lost its feasible incumbent')
    if progress_callback and upper < (current_score if workers == 1 else
                                      score_crossings(fixture, incumbent).crossings):
        progress_callback(optimized, upper)
    methods = {r['method'] for r in results}
    solver = ('exact-residual-factor-elimination' if methods <= {'residual-elimination'} else
              'monotone-component-branch-and-bound' if methods == {'implicit-branch-and-bound'} else
              'monotone-reduced-factor-pipeline')
    layout = ExactLayoutResult(initial, normalized, optimized,
        score_crossings(fixture, initial), score_crossings(fixture, normalized), score,
        sum(r['nodes'] for r in results), 'proven optimum' if lower == upper else 'bounded best known')
    diagnostics = [{**r['diagnostic'], 'method': r['method'], 'lower_bound': r['lower'],
                    'upper_bound': r['upper']} for r in results]
    return AutoLayoutResult(layout, solver, dict(
        lower_bound=lower, upper_bound=upper, optimality_gap=upper-lower,
        time_limit_seconds=time_limit_seconds, deadline_reached=expired(),
        component_workers=workers, component_diagnostics=diagnostics,
        prepared_reductions=[dict(component_id=i,
            chromosome_refs=sorted(r.label for r in bases[i].base_assignment),
            orientation_base={r.label: sign for r, sign in bases[i].base_assignment.items()},
            free_flip_groups=[[r.label for r in group] for group in bases[i].free_flip_groups],
            reason="disconnected incidence component and hard GF(2) equations",
            scope="component; shared by all stages") for i in range(len(components))],
        fallback_reason='; '.join(d.get('preflight_reason', '') for d in diagnostics),
        reduction_policy='no restart after factor reductions; reduced-factor continuation',
        warm_start='median homology-position sweep plus legal flips/local moves',
        conditional_bound_order_dp_max_chromosomes=12,
        heuristic_order_dp_max_chromosomes=12, heuristic_max_improving_steps=min(local_max_improving_steps, 20)))
