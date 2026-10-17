"""Mirror equivalence and selective neighborhood policy; independent controls."""
from pathlib import Path
import benchmark_efficiency as base

VARIANTS={
    'direct_control': ('milp_reclaim','weighted','highs'),
    'three_control': ('hybrid_3','weighted','highs'),
    'adaptive_control': ('hybrid_adaptive','weighted','highs'),
    'mirror_direct': ('milp_reclaim','weighted','highs',dict(mirror_symmetry=True)),
    'mirror_three': ('hybrid_3','weighted','highs',dict(mirror_symmetry=True)),
    'mirror_adaptive': ('hybrid_adaptive','weighted','highs',dict(mirror_symmetry=True)),
    'mirror_selective': ('budget_selective','weighted','highs',dict(mirror_symmetry=True,
        neighborhood_min_decisions=512,neighborhood_fraction=.15,
        neighborhood_max_seconds=30,deduplicate_neighborhoods=True)),
}


def design(entries):
    tasks=[]
    for budget,repeats,indices in ((150,1,range(9)),(600,2,range(6,9))):
        for repeat in range(repeats):
            for index in indices:
                for start in ('public','pre_hybrid'):
                    for variant in VARIANTS:
                        tasks.append(dict(case_index=index,case_id=entries[index]['case_id'],
                            budget=budget,repeat=repeat,start=start,variant=variant))
    return tasks


def prepare(args):
    if not args.reference_run:
        runs=[p for p in Path(args.history).glob('efficiency_*') if (p/'inputs.json').exists()
              and (p/'tasks.json').exists() and len(__import__('json').loads((p/'tasks.json').read_text()))==252]
        if not runs:raise ValueError('No reference efficiency run found; supply --reference-run')
        args.reference_run=str(max(runs,key=lambda p:(p.stat().st_mtime_ns,str(p))))
    original_prepare(args)

original_prepare=base.prepare

if __name__=='__main__':
    base.WORKER_ENTRY=Path(__file__).resolve()
    base.VARIANTS=VARIANTS
    base.design=design
    base.prepare=prepare
    base.main()
