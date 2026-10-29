"""Experimental pairwise interface solve followed by outward, one-sided refinement.

The entire propagated layout is one transaction against the full objective.
This is not a rigid-bundle contraction or a global optimality certificate.
"""
from dataclasses import replace
import math
import time
from .conditional_layout import PairScorer, row_candidates
from .heuristic import _component_map
from .layout import LayoutState, score_crossings, canonicalize_component_order
from .orientation_space import orientation_basis
from .saved_layout import validate_saved_layout, decode_saved_layout


def optimize_interface(fixture, starting_state, blocks=None, *, seconds=3600,
                       row_seconds=30, sweeps=1, progress=None, event=None,
                       boundary=None, pair_solver=None):
    if not math.isfinite(seconds) or seconds <= 0 or row_seconds <= 0:
        raise ValueError('Positive finite budgets required')
    if fixture.orientation_constraints:
        raise ValueError('Experimental interface propagation currently requires no hard orientation constraints')
    begin = time.monotonic()
    state = canonicalize_component_order(fixture, starting_state)
    basis = orientation_basis(fixture, fixture.chromosome_refs)
    validate_saved_layout(fixture, state, basis)
    species = list(fixture.species_ids)
    scorer = PairScorer(fixture)
    if boundary is None:
        index = max(range(len(scorer.pairs)), key=lambda i: score_crossings(scorer.fixtures[i], state).crossings)
    else:
        index = scorer.pairs.index(tuple(boundary))
    pair = scorer.pairs[index]
    sub = replace(fixture, chromosomes=tuple(c for c in fixture.chromosomes if c.ref.species_id in pair))
    seed = LayoutState({sp: state.chromosome_order[sp] for sp in pair},
                       {r: sign for r, sign in state.chromosome_orientation.items() if r.species_id in pair})
    if pair_solver is None:
        from .hybrid_experiment import run_hybrid
        pair_solver = run_hybrid
    if event: event(dict(stage='pair_solve', boundary=list(pair), elapsed=0))
    result = pair_solver(sub, seed, 'milp_reclaim', seconds * .65, backend='highs', scheduling='weighted')
    solved = decode_saved_layout(sub, result['optimized_state'])
    validate_saved_layout(sub, solved, orientation_basis(sub, sub.chromosome_refs))
    if score_crossings(sub, solved).crossings > score_crossings(sub, seed).crossings:
        raise ValueError('Pair solver regressed')
    order = dict(state.chromosome_order); order.update(solved.chromosome_order)
    signs = dict(state.chromosome_orientation); signs.update(solved.chromosome_orientation)
    proposal = LayoutState(order, signs)
    history = []
    # Fix each newly arranged row as the anchor for the next outward row.
    # All ambiguity stays in the canonical all-copy score; no chromosome mapping is invented.
    targets = [(j, j+1) for j in range(index-1, -1, -1)] + [(j, j-1) for j in range(index+2, len(species))]
    _, component_of = _component_map(fixture)
    deadline = begin + seconds
    for j, anchor in targets:
        if time.monotonic() >= deadline: break
        sp = species[j]; edge = min(j, anchor)
        before = score_crossings(scorer.fixtures[edge], proposal).crossings
        value = before; stop = min(deadline, time.monotonic()+row_seconds)
        if event: event(dict(stage='outward_propagation', species=sp, anchor=species[anchor], elapsed=time.monotonic()-begin))
        while time.monotonic() < stop:
            best = None; best_value = value
            for candidate in row_candidates(proposal, sp, component_of, basis):
                if time.monotonic() >= stop: break
                validate_saved_layout(fixture, candidate, basis)
                score = score_crossings(scorer.fixtures[edge], candidate).crossings
                if score < best_value: best, best_value = candidate, score
            if best is None: break
            proposal, value = best, best_value
        history.append(dict(species=sp, anchor=species[anchor], before=before, after=value))
    validate_saved_layout(fixture, proposal, basis)
    initial = score_crossings(fixture, state); candidate_score = score_crossings(fixture, proposal)
    accepted = candidate_score.crossings < initial.crossings
    final = proposal if accepted else state
    output = score_crossings(fixture, final)
    if output.link_count != initial.link_count: raise ValueError('Link count changed')
    if progress: progress(final)
    return dict(optimized_state=final.to_dict(), input=initial.to_dict(), output=output.to_dict(),
                status='interface-guided heuristic; no global certificate', global_lower_bound=None,
                boundary=list(pair), pair_result=result, propagation=history,
                proposal_state=proposal.to_dict(), proposal_score=candidate_score.to_dict(), accepted=accepted,
                seconds=time.monotonic()-begin,
                method='pairwise MILP followed by outward one-sided chromosome refinement; transactional full-score acceptance',
                limitation='Propagation is greedy, not a rigid bundle permutation; pair preparation may exceed requested solver budget')
