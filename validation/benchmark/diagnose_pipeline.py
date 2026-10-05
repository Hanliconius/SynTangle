"""Read saved audits or sample a time-limited solve without changing decisions."""
import argparse
from contextlib import ExitStack
import csv
import faulthandler
import json
from pathlib import Path
import signal
import time
from unittest.mock import patch


def atomic(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2))
    temporary.replace(path)


def audit(directory):
    roots = sorted(Path(directory).glob('pipeline_pair_*'))
    if not roots:
        raise ValueError(f'No pipeline_pair_* results under {directory}')
    for root in roots:
        for path in sorted(root.glob('*/paired.json')):
            payload = json.loads(path.read_text())
            print(f'CASE {payload["case_id"]} RUN {root.name}')
            for method in payload['methods']:
                print(f'  {method["method"]}: {method["status"]}')
                if method['method'] != 'pipeline':
                    continue
                result = method.get('result')
                if not result:
                    incumbent = method.get('incumbent', {})
                    print(f'  saved incumbent C={incumbent.get("crossings", "missing")} '
                          f'at {incumbent.get("seconds", "unknown")} seconds; '
                          'NO FINAL REDUCTION AUDIT OR PROOF GAP')
                    continue
                details = result['solver'] if 'component_diagnostics' not in result else result
                # AutoLayoutResult flattens its details into the serialized result.
                if not isinstance(details, dict):
                    raise ValueError(f'Missing solver diagnostics in {path}')
                print(f'  bounds={details.get("lower_bound")}..{details.get("upper_bound")} '
                      f'gap={details.get("optimality_gap")}')
                for component in details['component_diagnostics']:
                    print('  COMPONENT', component.get('component_id'),
                          'method=', component.get('method'),
                          'preflight=', component.get('preflight_reason', ''),
                          'nodes=', component.get('nodes_evaluated', component.get('continuation_nodes')),
                          'orientation_pruned=', component.get('orientation_branches_pruned'),
                          'orientation_forced=', component.get('orientation_groups_forced'))
                    for handoff in component.get('handoffs', []):
                        print('    HANDOFF', json.dumps(handoff, sort_keys=True))


class DiagnosticCutoff(Exception):
    pass


class Observer:
    """Aggregate call timing plus small proof/context samples, with exact returns."""
    def __init__(self):
        self.started = time.perf_counter()
        self.stats = {}
        self.stack = []
        self.samples = []
        self.incumbents = []
        self.patches = ExitStack()

    def sample(self, label, value):
        # Keep a bounded sample for each event type; counters retain full totals.
        if sum(s['event'] == label for s in self.samples) < 12:
            self.samples.append(dict(event=label, seconds=time.perf_counter()-self.started, **value))

    def wrap(self, module, name, label=None, before=None, after=None):
        original = getattr(module, name)
        label = label or name
        stats = self.stats.setdefault(label, dict(calls=0, completed=0, inclusive_seconds=0.,
                                                  exclusive_seconds=0., interrupted=0))
        def observed(*args, **kwargs):
            stats['calls'] += 1
            frame = dict(label=label, started=time.perf_counter(), children=0.)
            self.stack.append(frame)
            try:
                if before:
                    before(args, kwargs)
                result = original(*args, **kwargs)
                stats['completed'] += 1
                if after:
                    after(result, args, kwargs)
                return result
            except DiagnosticCutoff:
                stats['interrupted'] += 1
                raise
            except Exception as exc:
                self.sample(label+'.exception', dict(type=type(exc).__name__, message=str(exc)))
                raise
            finally:
                elapsed = time.perf_counter()-frame['started']
                stats['inclusive_seconds'] += elapsed
                stats['exclusive_seconds'] += max(0., elapsed-frame['children'])
                self.stack.pop()
                if self.stack:
                    self.stack[-1]['children'] += elapsed
        self.patches.enter_context(patch.object(module, name, observed))

    def install(self):
        import syntangle.search_pipeline as pipeline
        import syntangle.heuristic as heuristic
        import syntangle.branch_bound as branch
        import syntangle.residual_solver as residual
        def component(args, kwargs):
            self.sample('component.start', dict(component_id=args[3] if len(args)>10 else args[1]))
        self.wrap(pipeline, 'median_order_seed')
        self.wrap(heuristic, 'optimize_local_search')
        self.wrap(heuristic, 'optimize_species_component_order', label='heuristic.order_dp')
        self.wrap(pipeline, 'build_residual_factorization')
        self.wrap(pipeline, '_solve_incidence_component', before=component)
        self.wrap(pipeline, '_solve_branch_component', before=lambda a,k: self.sample(
            'implicit.start', dict(component_id=a[3], free_groups=len(a[-1].free_flip_groups))))
        self.wrap(residual, '_build_domains')
        self.wrap(residual, '_build_factor_tables')
        self.wrap(residual, '_eliminate_variable', after=lambda result,a,k: self.sample(
            'elimination.completed', dict(variable=a[0], remaining_scope=list(result[1].remaining_scope))))
        self.wrap(residual, '_drop_invariant_variables', after=lambda result,a,k: self.sample(
            'scope.reduction', dict(factor=a[0].factor_id, removed=result[1], remaining_scope=list(result[0].scope))))
        def continuation(args, kwargs):
            factors, domains, budget, _ = args
            active = sorted({v for f in factors for v in f.scope})
            retained = [r.variable_id for r in budget.reconstruction_stack]
            if set(active) & (set(retained) | set(budget.conditions)):
                raise AssertionError('Removed variable reintroduced in continuation')
            self.sample('continuation.start', dict(active_variables=active,
                retained_eliminations=retained, conditions=dict(budget.conditions),
                scope_reductions=budget.scope_reductions))
        self.wrap(residual, '_bounded_factor_search', before=continuation)
        self.wrap(branch, 'build_relaxed_crossing_bound')
        def reduction(result, args, kwargs):
            stats = self.stats['_reduce_orientation']
            stats['forced_total'] = stats.get('forced_total', 0) + result.forced
            stats['pruned_calls'] = stats.get('pruned_calls', 0) + int(result.pruned)
            self.sample('orientation.reduction', dict(forced=result.forced, pruned=result.pruned,
                unresolved=sum(x is None for x in result.bits), lower_bound=result.lower_bound))
        self.wrap(branch, '_reduce_orientation', after=reduction)
        self.wrap(branch, '_search_orders_for_orientation')
        self.wrap(branch, '_conditional_edge_minimum', after=lambda result,a,k: self.stats[
            '_conditional_edge_minimum'].__setitem__('cache_hits',
                self.stats['_conditional_edge_minimum'].get('cache_hits', 0)+int(result[2])))
        self.wrap(branch, 'build_pairwise_order_costs', label='bound.pairwise_costs')
        self.wrap(branch, 'solve_order_subset_dp', label='bound.order_dp')
        self.wrap(branch, '_edge_cost')
        self.wrap(branch._Budget, 'consume', label='branch.budget.consume', after=lambda result,a,k:
            self.stats['branch.budget.consume'].__setitem__('accepted',
                self.stats['branch.budget.consume'].get('accepted', 0)+int(result)))
        return self

    def close(self):
        self.patches.close()


def diagnose(args):
    from syntangle import load_validation_bundle, optimize_auto, score_crossings
    root = Path(args.root).resolve()
    with (root/'benchmark_manifest.tsv').open() as handle:
        entry = list(csv.DictReader(handle, delimiter='\t'))[args.index]
    output = Path(args.output)/entry['case_id']
    output.mkdir(parents=True, exist_ok=False)
    fixture = load_validation_bundle(root/entry['case_dir'])
    observer = Observer().install()
    status = 'failed'
    result_payload = None
    interrupted_stack = []
    def cutoff(signum, frame):
        interrupted_stack[:] = [f['label'] for f in observer.stack]
        raise DiagnosticCutoff()
    def checkpoint(state, crossings):
        if score_crossings(fixture, state).crossings != crossings:
            raise AssertionError('Incumbent fails canonical scoring')
        record = dict(crossings=crossings, seconds=time.perf_counter()-observer.started)
        observer.incumbents.append(record)
        atomic(output/'incumbent.json', {**record, 'state':state.to_dict(), 'status':'proof pending'})
        print('INCUMBENT', crossings, f'{record["seconds"]:.3f}s', flush=True)
    previous_handler = signal.signal(signal.SIGALRM, cutoff)
    faulthandler.dump_traceback_later(30, repeat=True)
    signal.setitimer(signal.ITIMER_REAL, args.seconds)
    print('DIAGNOSING', entry['case_id'], 'seconds=', args.seconds, flush=True)
    try:
        result = optimize_auto(fixture, transition_cap_per_component=250000,
            branch_node_cap_per_component=1000, local_restarts=1, local_max_improving_steps=5,
            component_workers=1, seed=int(entry['seed']), progress_callback=checkpoint)
        result_payload = result.to_dict()
        status = 'complete'
    except DiagnosticCutoff:
        status = 'diagnostic cutoff; solve unresolved'
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)
        faulthandler.cancel_dump_traceback_later()
        observer.close()
        payload = dict(case_id=entry['case_id'], status=status,
            seconds=time.perf_counter()-observer.started, stats=observer.stats,
            interrupted_stack=interrupted_stack, samples=observer.samples,
            incumbents=observer.incumbents, result=result_payload,
            note='Diagnostic instrumentation adds overhead; these are not benchmark timings.')
        atomic(output/'diagnostic.json', payload)
        print('STATUS', status, flush=True)
        print('ACTIVE AT CUTOFF', ' -> '.join(interrupted_stack), flush=True)
        print('TIMING (exclusive seconds exclude other instrumented calls)', flush=True)
        for name, stats in sorted(observer.stats.items(), key=lambda item: -item[1]['exclusive_seconds']):
            print(name, json.dumps(stats, sort_keys=True), flush=True)
        print('REDUCTION SAMPLES', json.dumps(observer.samples, sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--audit')
    parser.add_argument('--root')
    parser.add_argument('--output')
    parser.add_argument('--index', type=int)
    parser.add_argument('--seconds', type=int, default=90)
    args = parser.parse_args()
    audit(args.audit) if args.audit else diagnose(args)


if __name__ == '__main__':
    main()
