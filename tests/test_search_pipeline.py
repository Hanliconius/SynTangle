"""Exhaustive correctness checks and adversarial solver-transition tests."""
import itertools
from pathlib import Path
import random
import unittest
from unittest.mock import patch

from syntangle import load_fixture, optimize_auto, score_crossings
from syntangle.residual_solver import (
    _Domain, _TableFactor, _WorkBudget, _MutableDiagnostics,
    _solve_factor_system, _drop_invariant_variables)

ROOT = Path(__file__).resolve().parents[1]


def solve(factors, domains, work_cap, node_cap, separator=8):
    budget = _WorkBudget(work_cap, 10000, continuation_node_cap=node_cap)
    diagnostics = _MutableDiagnostics()
    result = _solve_factor_system(factors, domains, budget, diagnostics,
                                 separator_domain_cap=separator)
    return result, budget, diagnostics


def exhaustive(factors, domains):
    variables = sorted({v for f in factors for v in f.scope})
    def score(assignment):
        return sum(f.values[tuple(assignment[v] for v in f.scope)] for f in factors)
    optimum = min(score(dict(zip(variables, values))) for values in
                  itertools.product(*(range(domains[v].size) for v in variables)))
    return optimum, score


class PipelineTests(unittest.TestCase):
    def test_bounds_and_reconstruction_under_random_budget_interruptions(self):
        rng = random.Random(45)
        for case in range(30):
            domains = {v: _Domain(v, (0, 1)) for v in 'abcde'}
            scopes = [('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'e'), ('e', 'b')]
            factors = [_TableFactor(str(i), scope,
                       {key: rng.randrange(10) for key in itertools.product((0, 1), repeat=2)})
                       for i, scope in enumerate(scopes)]
            optimum, score = exhaustive(factors, domains)
            for work in (0, 4, 8, 10000):
                for nodes in (0, 1, 1000):
                    result, budget, _ = solve(factors, domains, work, nodes)
                    self.assertEqual(score(result.assignment), result.cost)
                    self.assertLessEqual(result.lower, optimum)
                    self.assertGreaterEqual(result.cost, optimum)
                    if nodes == 1000:
                        self.assertEqual(result.lower, optimum)
                        self.assertEqual(result.cost, optimum)
                    for handoff in budget.handoffs:
                        active = set(handoff['active_variables'])
                        self.assertFalse(active & set(handoff['retained_conditions']))
                        self.assertFalse(active & {r['variable'] for r in handoff['retained_eliminations']})

    def test_eliminated_leaf_not_reintroduced(self):
        domains = {v: _Domain(v, (0, 1)) for v in 'abc'}
        factors = [_TableFactor('ab', ('a', 'b'), {(0,0): 2, (0,1): 0, (1,0): 0, (1,1): 2}),
                   _TableFactor('bc', ('b', 'c'), {(0,0): 3, (0,1): 0, (1,0): 0, (1,1): 3})]
        result, budget, _ = solve(factors, domains, 4, 100)
        self.assertTrue(budget.handoffs)
        transition = budget.handoffs[0]
        self.assertEqual([r['variable'] for r in transition['retained_eliminations']], ['a'])
        self.assertNotIn('a', transition['active_variables'])
        self.assertEqual(result.cost, 0)
        self.assertEqual(result.lower, 0)
        self.assertIn('a', result.assignment)

    def test_conditional_reductions_remain_scoped_to_branch(self):
        domains = {v: _Domain(v, (0, 1)) for v in 'abcde'}
        scopes = [('a','b'), ('b','c'), ('c','a'), ('c','d'), ('d','e'), ('e','c')]
        factors = [_TableFactor(str(i), scope, {key: int(key[0] == key[1])
                  for key in itertools.product((0,1), repeat=2)}) for i, scope in enumerate(scopes)]
        result, budget, _ = solve(factors, domains, 0, 100)
        conditions = [h['retained_conditions'] for h in budget.handoffs]
        self.assertTrue(any(c.get('c') == 0 for c in conditions))
        self.assertTrue(any(c.get('c') == 1 for c in conditions))
        optimum, score = exhaustive(factors, domains)
        self.assertEqual(result.cost, optimum)
        self.assertEqual(result.lower, optimum)
        self.assertEqual(score(result.assignment), optimum)

    def test_scope_invariance_and_disconnected_piece_preserved(self):
        domains = {v: _Domain(v, (0, 1)) for v in 'abc'}
        factor = _TableFactor('invariant', ('a','b'),
                              {(0,0): 3, (1,0): 3, (0,1): 1, (1,1): 1})
        reduced, removed = _drop_invariant_variables(factor, domains)
        self.assertEqual(removed, 1)
        self.assertEqual(reduced.scope, ('b',))
        result, budget, _ = solve([reduced, _TableFactor('c', ('c',), {(0,):2, (1,):4})],
                                  domains, 2, 10)
        self.assertEqual(result.cost, 3)
        self.assertTrue(budget.handoffs)
        self.assertNotIn('a', budget.handoffs[0]['active_variables'])
        self.assertEqual(len(budget.handoffs[0]['active_variables']), 1)
        self.assertEqual(budget.handoffs[0]['leaf_eliminations_retained'], 1)

    def test_preflight_fails_before_any_table_scope_reduction(self):
        from syntangle.heuristic import _component_map
        from syntangle.residual import build_residual_factorization
        from syntangle.residual_solver import _solve_incidence_component
        from syntangle.layout import SearchSpaceTooLarge
        fixture = load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        components, component_of = _component_map(fixture)
        graph = build_residual_factorization(fixture)
        with patch('syntangle.residual_solver._drop_invariant_variables',
                   side_effect=AssertionError('reduction before preflight')):
            with self.assertRaises(SearchSpaceTooLarge):
                _solve_incidence_component(fixture, 0, components[0], component_of,
                    tuple(v for v in graph.variables if v.component_id == 0),
                    tuple(f for f in graph.factors if f.component_id == 0),
                    40320, 1, 1, 8, continuation_node_cap=10)

    def test_exhausted_node_budget_does_not_claim_proof(self):
        domains = {v: _Domain(v, (0,1)) for v in 'abc'}
        factors = [_TableFactor(str(i), scope, {key:int(key[0] == key[1])
                  for key in itertools.product((0,1), repeat=2)})
                  for i, scope in enumerate([('a','b'), ('b','c'), ('c','a')])]
        result, budget, _ = solve(factors, domains, 0, 0)
        self.assertEqual(result.cost, 3)
        self.assertEqual(result.lower, 0)
        self.assertTrue(budget.handoffs[0]['unresolved'])
        optimum, _ = exhaustive(factors, domains)
        self.assertEqual(optimum, 1)

    def test_completed_component_not_repeated_when_later_component_dispatches(self):
        from syntangle.model import Fixture, Chromosome, ChromosomeRef, BlockOccurrence
        import syntangle.search_pipeline as pipeline
        chromosomes = []
        for species in ('A', 'B'):
            chromosomes.append(Chromosome(ChromosomeRef(species, 'solo'), 10, 0, 1,
                (BlockOccurrence(species+'solo', 'solo', 1, 2, '+'),)))
            for i in range(3):
                blocks = tuple(BlockOccurrence(f'{species}{i}_{j}',
                    f'large_{i}_{j}' if species == 'A' else f'large_{j}_{i}',
                    j*2+1, j*2+2, '+') for j in range(3))
                chromosomes.append(Chromosome(ChromosomeRef(species, str(i)), 10, i+1, 1, blocks))
        fixture = Fixture(1, 'mixed', 'mixed', 'mixed dispatch', tuple(chromosomes), (), {})
        with patch.object(pipeline, '_solve_incidence_component', wraps=pipeline._solve_incidence_component) as exact:
            with patch.object(pipeline, '_solve_branch_component', wraps=pipeline._solve_branch_component) as bounded:
                result = optimize_auto(fixture, permutation_cap_per_species=2, local_restarts=1,
                                       branch_node_cap_per_component=1000)
        self.assertEqual(exact.call_count, 2)
        self.assertEqual(bounded.call_count, 1)
        self.assertEqual(len(result.details['component_diagnostics']), 2)
        # The identical propagated basis object reaches both preparation and dispatch.
        failed = next(call for call in exact.call_args_list if len(call.args[4]) > 2)
        self.assertIs(failed.kwargs['prepared_basis'], bounded.call_args.args[-1])

    def test_auto_does_not_restart_legacy_exact_methods(self):
        fixture = load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        with patch('syntangle.layer_dp.exact_optimize_layer_dp', side_effect=AssertionError('restarted')):
            result = optimize_auto(fixture, transition_cap_per_component=1,
                                   branch_node_cap_per_component=1000, local_restarts=1)
        self.assertEqual(score_crossings(fixture, result.layout.optimized_state).crossings,
                         result.details['upper_bound'])
        self.assertTrue(all(d['factor_reductions_before_dispatch'] == 0
                            for d in result.details['component_diagnostics']))

    def test_large_order_relaxation_is_admissible_without_subset_dp(self):
        from syntangle.branch_bound import _conditional_edge_minimum, _edge_cost
        from syntangle.heuristic import _component_map
        from syntangle.layout import initial_layout_state, LayoutState
        fixture = load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        components, _ = _component_map(fixture)
        state = initial_layout_state(fixture)
        checks = 0
        for nodes in components:
            for target_index in range(len(fixture.species_ids)):
                for neighbor_index in (target_index-1, target_index+1):
                    if not 0 <= neighbor_index < len(fixture.species_ids):
                        continue
                    with patch('syntangle.branch_bound.solve_order_subset_dp',
                               side_effect=AssertionError('subset DP should be skipped')):
                        bound, preferred, _ = _conditional_edge_minimum(fixture, state,
                            target_index, neighbor_index, nodes, cache={}, order_dp_max_chromosomes=0)
                    target, neighbor = fixture.species_ids[target_index], fixture.species_ids[neighbor_index]
                    from syntangle.incidence import chromosome_node_id
                    other = tuple(r for r in state.chromosome_order[neighbor] if chromosome_node_id(r) in nodes)
                    for permutation in itertools.permutations(preferred):
                        orders = dict(state.chromosome_order)
                        orders[target] = permutation
                        candidate = LayoutState(orders, state.chromosome_orientation)
                        if target_index < neighbor_index:
                            actual = _edge_cost(fixture, target, neighbor, permutation, other, candidate, nodes)
                        else:
                            actual = _edge_cost(fixture, neighbor, target, other, permutation, candidate, nodes)
                        self.assertLessEqual(bound, actual)
                        checks += 1
        self.assertGreater(checks, 0)

    def test_actual_layout_continues_after_leaf_eliminations(self):
        from syntangle.model import Fixture, Chromosome, ChromosomeRef, BlockOccurrence
        chromosomes = tuple(Chromosome(ChromosomeRef(f'sp{i}', 'chr'), 10, 0, 1,
            tuple(BlockOccurrence(f'{i}_{j}', f'H{j}', j*2+1, j*2+2, '+') for j in range(3)))
            for i in range(8))
        fixture = Fixture(1, 'chain', 'chain', 'continuation', chromosomes, (), {})
        result = optimize_auto(fixture, transition_cap_per_component=40,
                               branch_node_cap_per_component=10, local_restarts=1)
        self.assertEqual(result.solver, 'monotone-reduced-factor-pipeline')
        handoffs = result.details['component_diagnostics'][0]['handoffs']
        self.assertTrue(handoffs)
        self.assertTrue(handoffs[0]['retained_eliminations'])
        self.assertEqual(result.details['upper_bound'], 0)
        self.assertEqual(score_crossings(fixture, result.layout.optimized_state).crossings, 0)
        self.assertEqual(result.layout.optimality_status, 'proven optimum')

    def test_real_fixture_continuation_bound_matches_unlimited_optimum(self):
        fixture = load_fixture(ROOT/'examples/fixtures/fusion_chain_closed_cycle.json')
        exact = optimize_auto(fixture, local_restarts=1)
        optimum = exact.layout.optimized_score.crossings
        for cap in (20, 40, 80, 160):
            result = optimize_auto(fixture, transition_cap_per_component=cap,
                                   branch_node_cap_per_component=1000, local_restarts=1)
            self.assertLessEqual(result.details['lower_bound'], optimum)
            self.assertGreaterEqual(result.details['upper_bound'], optimum)
            self.assertEqual(result.details['upper_bound'], optimum)


if __name__ == '__main__':
    unittest.main()
