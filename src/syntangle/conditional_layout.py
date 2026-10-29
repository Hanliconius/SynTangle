"""Bounded conditional layout improvement without constructing a global MILP.

Temporary fixed rows are neighborhood definitions, never proof-based exclusions.
Every crossing and chromosome remains in the objective and reconstructed state.
"""
from dataclasses import replace
from itertools import product
import math
import time
from .layout import LayoutState, score_crossings, canonicalize_component_order
from .heuristic import _component_map
from .orientation_space import orientation_basis
from .saved_layout import validate_saved_layout


class PairScorer:
    """Canonical inversion scorer on only the original adjacent pairs touched."""
    def __init__(self, fixture):
        self.fixture = fixture
        self.pairs = tuple(zip(fixture.species_ids, fixture.species_ids[1:]))
        self.fixtures = [replace(fixture, chromosomes=tuple(c for c in fixture.chromosomes
                         if c.ref.species_id in pair), orientation_constraints=())
                         for pair in self.pairs]

    def score(self, state, active=None):
        return sum(score_crossings(f, state).crossings
                   for pair, f in zip(self.pairs, self.fixtures)
                   if active is None or set(pair) & set(active))


def mirror_block(state, species, component_of):
    order = dict(state.chromosome_order)
    signs = dict(state.chromosome_orientation)
    for sp in species:
        grouped = {}
        for ref in order[sp]:
            grouped.setdefault(component_of[ref], []).append(ref)
            signs[ref] *= -1
        order[sp] = tuple(ref for group in grouped.values() for ref in reversed(group))
    return LayoutState(order, signs)


def align_blocks(fixture, state, blocks, scorer, component_of):
    """Exact binary DP over these block mirrors only; not global optimality."""
    if [sp for block in blocks for sp in block] != list(fixture.species_ids):
        raise ValueError('Blocks must partition the displayed species in order')
    if fixture.orientation_constraints:
        return state, {'status': 'mirror stage skipped: hard constraints require joint handling'}
    alternatives = [(state, mirror_block(state, block, component_of)) for block in blocks]
    unary = []
    for block, states in zip(blocks, alternatives):
        internal = [i for i, pair in enumerate(scorer.pairs) if set(pair) <= set(block)]
        unary.append([sum(score_crossings(scorer.fixtures[i], st).crossings for i in internal)
                      for st in states])
    dp = [(unary[0][bit], (bit,)) for bit in (0, 1)]
    for i in range(1, len(blocks)):
        boundary = scorer.pairs.index((blocks[i-1][-1], blocks[i][0]))
        choices = []
        for bit in (0, 1):
            options = []
            for prev in (0, 1):
                st = state
                if prev: st = mirror_block(st, blocks[i-1], component_of)
                if bit: st = mirror_block(st, blocks[i], component_of)
                cost = score_crossings(scorer.fixtures[boundary], st).crossings
                options.append((dp[prev][0] + unary[i][bit] + cost, dp[prev][1] + (bit,)))
            choices.append(min(options))
        dp = choices
    value, bits = min(dp)
    result = state
    for block, bit in zip(blocks, bits):
        if bit: result = mirror_block(result, block, component_of)
    if scorer.score(result) != value:
        raise ValueError('Block DP disagrees with full crossing objective')
    return result, {'status': 'exact optimum within two-mirror-per-block domain only',
                    'mirror_bits': list(bits), 'crossings': value}


def row_candidates(state, species, component_of, basis):
    """Complete-chromosome relocations, reversals and legal orientation groups."""
    refs = state.chromosome_order[species]
    for i, ref in enumerate(refs):
        positions = [j for j, other in enumerate(refs) if component_of[other] == component_of[ref]]
        for j in positions:
            if i == j: continue
            row = list(refs); row.insert(j, row.pop(i))
            order = dict(state.chromosome_order); order[species] = tuple(row)
            yield LayoutState(order, dict(state.chromosome_orientation))
    for group in basis.free_flip_groups:
        if all(ref.species_id == species for ref in group):
            signs = dict(state.chromosome_orientation)
            for ref in group: signs[ref] *= -1
            yield LayoutState(dict(state.chromosome_order), signs)
    # Joint order+orientation reversal can escape single-move local minima.
    if not any(ref.species_id != species for group in basis.free_flip_groups
               if any(ref.species_id == species for ref in group) for ref in group):
        yield mirror_block(state, [species], component_of)


def optimize_conditional(fixture, starting_state, blocks, *, seconds=3600,
                         row_seconds=20, sweeps=3, progress=None, event=None):
    """Align genus blocks, then boundary and overlapping row neighborhoods.

    Fixed neighbors contribute their exact costs. Strict improvement only.
    This is a time-bounded conditional heuristic, with no global lower bound.
    """
    if not math.isfinite(seconds) or seconds <= 0 or not math.isfinite(row_seconds) or row_seconds <= 0 or sweeps < 1:
        raise ValueError('Positive finite budgets and sweeps required')
    begin = time.monotonic(); deadline = begin + seconds
    basis = orientation_basis(fixture, fixture.chromosome_refs)
    validate_saved_layout(fixture, starting_state, basis)
    state = canonicalize_component_order(fixture, starting_state)
    _, component_of = _component_map(fixture)
    scorer = PairScorer(fixture); initial = score_crossings(fixture, starting_state)
    if event: event({'stage': 'block_alignment', 'elapsed': time.monotonic()-begin})
    state, block_info = align_blocks(fixture, state, blocks, scorer, component_of)
    if scorer.score(state) > initial.crossings:
        raise ValueError('Block alignment lost incumbent')
    if progress: progress(state)
    boundary_rows = list(dict.fromkeys(sp for left, right in zip(blocks, blocks[1:])
                                      for sp in (left[-1], right[0])))
    history = []; evaluations = 0
    for sweep in range(sweeps):
        # First focus on genus boundaries; subsequent sweeps release every row.
        rows = boundary_rows if sweep == 0 else list(fixture.species_ids)
        if sweep % 2: rows.reverse()
        for sp in rows:
            if time.monotonic() >= deadline: break
            stop = min(deadline, time.monotonic() + row_seconds)
            before = scorer.score(state, [sp]); value = before
            if event: event({'stage': 'row_neighborhood', 'sweep': sweep, 'species': sp,
                             'elapsed': time.monotonic()-begin})
            while time.monotonic() < stop:
                best = None; best_value = value
                for candidate in row_candidates(state, sp, component_of, basis):
                    if time.monotonic() >= stop: break
                    try: validate_saved_layout(fixture, candidate, basis)
                    except ValueError: continue
                    score = scorer.score(candidate, [sp]); evaluations += 1
                    if score < best_value: best, best_value = candidate, score
                if best is None: break
                state, value = best, best_value
                if progress: progress(state)
            history.append({'sweep': sweep, 'species': sp, 'before_affected': before,
                            'after_affected': value, 'outside_rows_fixed': True})
        # Overlapping multi-row mirrors provide joint moves after boundary sweeps.
        if sweep > 0:
            width = 3 if sweep == 1 else 5
            for left in range(max(1, len(fixture.species_ids)-width+1)):
                if time.monotonic() >= deadline: break
                active = list(fixture.species_ids[left:left+width])
                candidate = mirror_block(state, active, component_of)
                try: validate_saved_layout(fixture, candidate, basis)
                except ValueError: continue
                before = scorer.score(state, active); after = scorer.score(candidate, active)
                evaluations += 1
                if after < before:
                    state = candidate
                    if progress: progress(state)
                history.append({'sweep': sweep, 'joint_window': active,
                                'before_affected': before, 'after_affected': min(before,after),
                                'outside_rows_fixed': True})
        if time.monotonic() >= deadline: break
    validate_saved_layout(fixture, state, basis)
    final = score_crossings(fixture, state)
    if final.crossings > initial.crossings or final.link_count != initial.link_count:
        raise ValueError('Conditional refinement regressed or changed links')
    return {'optimized_state': state.to_dict(), 'input': initial.to_dict(),
            'output': final.to_dict(), 'seconds': time.monotonic()-begin,
            'status': 'conditional heuristic; no global optimality certificate',
            'global_lower_bound': None, 'block_alignment': block_info,
            'row_history': history, 'candidate_evaluations': evaluations,
            'fixed_scope': 'temporary neighborhood only; never exported as forced decisions',
            'method': 'binary block-mirror DP then exact-scored whole-chromosome row moves'}
