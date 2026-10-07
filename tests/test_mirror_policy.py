from itertools import product,permutations
from pathlib import Path
import unittest
import numpy as np
from test_coupled_bound import fixture_for
from syntangle.global_experiment import JointModel
from syntangle.hybrid_experiment import run_hybrid,improve_neighborhoods
from syntangle.layout import LayoutState,initial_layout_state,score_crossings
from syntangle.incidence import chromosome_node_id
from syntangle.model import Fixture,Chromosome,ChromosomeRef,BlockOccurrence


def model_for(f):
    return JointModel(f,frozenset(chromosome_node_id(r) for r in f.chromosome_refs),f.chromosome_refs)

class MirrorPolicyTests(unittest.TestCase):
    def test_exact_mirror_equivalence_and_both_anchor_optima(self):
        f=fixture_for(19);m=model_for(f)
        raw=initial_layout_state(f);start=LayoutState(raw.chromosome_order,m.basis.base_assignment)
        self.assertTrue(m.mirror_certificate()['verified'])
        optima=[float('inf'),float('inf')]
        for bits in product((0,1),repeat=m.n):
            bits=np.asarray(bits);state=m.decode(bits,start)
            reflected=m.decode(1-bits,state)
            self.assertEqual(m.objective(bits),m.objective(1-bits))
            self.assertEqual(score_crossings(f,state).crossings,score_crossings(f,reflected).crossings)
            self.assertEqual(state,m.decode(1-m.encode(reflected),reflected))
            optima[bits[0]]=min(optima[bits[0]],m.objective(bits))
        self.assertEqual(optima[0],optima[1])
        for candidate in (start,m.decode(1-m.encode(start),start)):
            state,bound,info=m.highs(candidate,5,mirror_symmetry=True)
            self.assertEqual(bound,optima[0])
            anchor=info['mirror_reduction']['anchor_index']
            self.assertEqual(m.encode(state)[anchor],m.encode(candidate)[anchor])
            self.assertTrue(info['mirror_reduction']['applied'])
            # Independent audit keeps the original unanchored model.
            self.assertTrue(m.highs(state,5,strict_proof=True)['no_better_proven'])

    def test_three_chromosome_total_order_complement_remains_legal(self):
        chromosomes=[]
        for sp,groups in [('A',((0,1),(2,3),(4,5))),('B',((1,3),(5,0),(2,4)))]:
            for i,group in enumerate(groups):
                blocks=tuple(BlockOccurrence(f'{sp}_{h}',f'H{h}',k*3,k*3+1,'+') for k,h in enumerate(group))
                chromosomes.append(Chromosome(ChromosomeRef(sp,str(i)),20,i,1,blocks))
        f=Fixture(1,'mirror3','mirror3','three chromosome mirror',tuple(chromosomes),(),{})
        m=model_for(f);start=initial_layout_state(f)
        self.assertTrue(m.mirror_certificate()['verified'])
        self.assertTrue(m.triangles)
        for orders in product(*(permutations(row) for row in m.rows.values())):
            state=LayoutState(dict(zip(f.species_ids,orders)),m.basis.base_assignment)
            bits=m.encode(state);reflected=m.decode(1-bits,state)
            self.assertEqual(score_crossings(f,state).crossings,score_crossings(f,reflected).crossings)
        self.assertEqual(m.highs(start,5)[1],m.highs(start,5,mirror_symmetry=True)[1])

    def test_asymmetric_and_noninteger_models_do_not_authorize_pruning(self):
        m=model_for(fixture_for(19));m.linear[0]+=1
        self.assertFalse(m.mirror_certificate()['verified'])
        m.linear[0]-=.5
        self.assertFalse(m.mirror_certificate()['verified'])

    def test_selective_threshold_skips_heuristic_and_keeps_exact_optimum(self):
        f=fixture_for(31);m=model_for(f)
        raw=initial_layout_state(f);start=LayoutState(raw.chromosome_order,m.basis.base_assignment)
        r=run_hybrid(f,start,'hybrid_adaptive',5,mirror_symmetry=True,
            neighborhood_min_decisions=512,neighborhood_fraction=.15,neighborhood_max_seconds=30,
            deduplicate_neighborhoods=True)
        self.assertEqual(r['upper_bound'],r['lower_bound'])
        self.assertEqual(r['global_search_restarts'],0)
        self.assertTrue(all('neighborhood_search' not in d for d in r['component_diagnostics']))
        self.assertTrue(any(d.get('mirror_reduction',{}).get('applied') for d in r['component_diagnostics']))

    def test_duplicate_edge_windows_are_skipped_without_freezing_global_choices(self):
        f=fixture_for(19);m=model_for(f)
        raw=initial_layout_state(f);start=LayoutState(raw.chromosome_order,m.basis.base_assignment)
        state,info=improve_neighborhoods(m,start,2,deduplicate=True)
        self.assertLessEqual(info['subsolves'],2)
        self.assertTrue(info['unique_windows_per_sweep'])
        self.assertEqual(m.highs(state,5)[1],m.highs(start,5)[1])

class PolicyDesignTests(unittest.TestCase):
    def test_balanced_controls_and_repeated_long_runs(self):
        import sys
        sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'validation/benchmark'))
        import benchmark_policy as p
        tasks=p.design([dict(case_id=f'case_{i}') for i in range(9)])
        self.assertEqual(len(tasks),210)
        self.assertEqual(len({(r['case_id'],r['start'],r['variant'],r['budget'],r['repeat']) for r in tasks}),210)
        self.assertEqual(sum(r['budget']==600 for r in tasks),84)
        self.assertEqual(sum(r['budget']==150 for r in tasks),126)

    def test_policy_task_launches_policy_worker_and_collector(self):
        import json,os,subprocess,sys,tempfile
        from syntangle import load_validation_bundle
        import benchmark_efficiency as runner
        repo=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);case=root/'tiny';case.mkdir();out=root/'output';out.mkdir()
            (case/'species.tsv').write_text('species_id\tspecies_rank\nA\t1\nB\t2\n')
            (case/'chromosomes.tsv').write_text('species_id\tchromosome_id\tlength\tdisplay_rank\nA\t1\t10\t1\nB\t1\t10\t1\n')
            (case/'occurrences.tsv').write_text('occurrence_id\thomology_id\tspecies_id\tchromosome_id\tstart\tend\tstrand\nA1\tg1\tA\t1\t0\t1\t+\nA2\tg2\tA\t1\t3\t4\t+\nB2\tg2\tB\t1\t0\t1\t+\nB1\tg1\tB\t1\t3\t4\t+\n')
            (root/'benchmark_manifest.tsv').write_text('case_id\tcase_dir\tseed\ntiny\ttiny\t1\n')
            fixture=load_validation_bundle(case);state=initial_layout_state(fixture)
            (out/'tiny.public.json').write_text(json.dumps(dict(state=state.to_dict(),score=score_crossings(fixture,state).crossings)))
            (out/'inputs.json').write_text(json.dumps([dict(case_id='tiny',input_files=runner.fingerprint(case))]))
            (out/'tasks.json').write_text(json.dumps([dict(case_index=0,case_id='tiny',start='public',variant='mirror_selective',budget=2,repeat=0)]))
            env={**os.environ,'PYTHONPATH':str(repo/'src')}
            command=[sys.executable,str(repo/'validation/benchmark/benchmark_policy.py')]
            run=subprocess.run(command+['--root',str(root),'--output',str(out),'--index','0'],env=env,capture_output=True,text=True,timeout=15)
            self.assertEqual(run.returncode,0,run.stderr)
            saved=out/'tiny'/'mirror_selective_public_2_r0'/'task.json'
            result=json.loads(saved.read_text())['result']
            self.assertEqual(result['upper_bound'],result['lower_bound'])
            self.assertEqual(result['profile']['mirror_reduced_components'],1)
            collect=subprocess.run(command+['--collect','--output',str(out)],env=env,capture_output=True,text=True,timeout=10)
            self.assertEqual(collect.returncode,0,collect.stderr)
            self.assertIn('Expected 1 tasks; found 1',collect.stdout)
