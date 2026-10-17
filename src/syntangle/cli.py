from __future__ import annotations

import argparse
import json
from pathlib import Path

from .audit import build_layout_audit
from .branch_bound import optimize_branch_and_bound
from .decomposition import decompose_incidence_graph
from .fixtures import load_fixture
from .heuristic import optimize_auto, optimize_local_search
from .incidence import analyze_fixture, build_incidence_graph
from .layer_dp import exact_optimize_layer_dp
from .layout import exact_optimize_small
from .ordering import derive_ordering_constraints
from .orientation import solve_orientation_constraints
from .structural import build_structural_projection
from .visualize import write_layout_comparison_svg


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="syntangle")
    from importlib.metadata import version, PackageNotFoundError
    try: package_version = version('syntangle')
    except PackageNotFoundError: package_version = 'source-checkout'
    parser.add_argument('--version', action='version', version=package_version)
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser(
        "analyze-fixture",
        help="parse a synthetic fixture and summarize its constrained structure",
    )
    analyze.add_argument("fixture")
    analyze.add_argument("--pretty", action="store_true")

    optimize = subparsers.add_parser(
        "optimize-fixture",
        help="minimize legal crossing tangledness",
    )
    optimize.add_argument("fixture")
    optimize.add_argument(
        "--solver",
        choices=("auto", "branch-bound", "layer-dp", "enumerate", "local-search"),
        default="auto",
    )
    optimize.add_argument("--state-cap-per-component", type=int, default=250000)
    optimize.add_argument("--orientation-cap-per-component", type=int, default=4096)
    optimize.add_argument("--permutation-cap-per-species", type=int, default=40320)
    optimize.add_argument("--transition-cap-per-component", type=int, default=5000000)
    optimize.add_argument("--branch-node-cap-per-component", type=int, default=250000)
    optimize.add_argument("--local-restarts", type=int, default=8)
    optimize.add_argument("--local-max-improving-steps", type=int, default=10000)
    optimize.add_argument("--seed", type=int, default=1)
    optimize.add_argument("--audit-json")
    optimize.add_argument("--svg")
    optimize.add_argument("--pretty", action="store_true")
    refine = subparsers.add_parser('refine-genespace', help='Refine a complete exported layout on fixed block evidence')
    refine.add_argument('fixture')
    refine.add_argument('--layout', required=True, help='TSV: species, chromosome, rank, orientation')
    refine.add_argument('--output', required=True)
    refine.add_argument('--mode', choices=('flips','full'), default='flips')
    refine.add_argument('--seconds', type=float, default=60)
    ribbon = subparsers.add_parser('refine-ribbons', help='Refine whole-chromosome layout with explicit ribbon-twist objective')
    ribbon.add_argument('fixture')
    ribbon.add_argument('--layout', required=True, help='Complete TSV: species, chromosome, rank, orientation')
    ribbon.add_argument('--output', required=True)
    ribbon.add_argument('--objective', choices=('crossing-cap', 'combined'), default='crossing-cap')
    ribbon.add_argument('--twist-metric', choices=('count', 'width'), default='width')
    ribbon.add_argument('--weight', type=float, default=1.0)
    ribbon.add_argument('--seconds', type=float, default=600)
    for command, help_text in [('flip-assist','Whole-chromosome flips with fixed order'), ('run','Full chromosome ordering and orientation optimization')]:
        public = subparsers.add_parser(command, help=help_text)
        public.add_argument('--fixture')
        public.add_argument('--chromosomes')
        public.add_argument('--links')
        public.add_argument('--layout')
        public.add_argument('--output',required=True)
        public.add_argument('--seconds',type=float,default=60 if command=='flip-assist' else 600)
        public.add_argument('--figure',action='store_true')
        public.add_argument('--seed',type=int,default=1)
        if command=='run':
            public.add_argument('--ribbon-objective',choices=('none','crossing-cap','combined'),default='none')
            public.add_argument('--twist-metric',choices=('width','count'),default='width')
            public.add_argument('--weight',type=float,default=1)
            public.add_argument('--ribbon-seconds',type=float,default=600)
    for command in ('import-paf','import-jcvi'):
        importer = subparsers.add_parser(command, help='Experimental external-format conversion; not validated on real datasets')
        for name in ('query-species','target-species','query-sizes','target-sizes','output'):
            importer.add_argument('--'+name,required=True)
        if command=='import-paf':
            importer.add_argument('--paf',required=True)
            importer.add_argument('--min-mapq',type=int,default=0)
            importer.add_argument('--min-alignment-length',type=int,default=0)
        else:
            importer.add_argument('--simple',required=True,help='JCVI six-column .simple blocks')
            importer.add_argument('--query-bed',required=True)
            importer.add_argument('--target-bed',required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command in ('import-paf','import-jcvi'):
        from .importers import run_import
        return run_import(args)

    if args.command in ('flip-assist','run'):
        from .workflow import run
        return run(args)

    if args.command == 'refine-ribbons':
        from .ribbon_objective import refine_ribbons
        return refine_ribbons(args)

    if args.command == 'refine-genespace':
        if args.seconds <= 0:raise ValueError('--seconds must be positive')
        from .genespace import refine
        return refine(args)

    if args.command == "analyze-fixture":
        fixture = load_fixture(args.fixture)
        graph = build_incidence_graph(fixture)
        output = analyze_fixture(fixture).to_dict()
        output["orientation"] = solve_orientation_constraints(fixture).to_dict()
        output["decomposition"] = decompose_incidence_graph(graph).to_dict()
        output["ordering_constraints"] = derive_ordering_constraints(fixture).to_dict()
        output["structural_projection"] = build_structural_projection(fixture).to_dict()
        print(json.dumps(output, indent=2 if args.pretty else None, sort_keys=True))
        return 0

    fixture = load_fixture(args.fixture)

    if args.solver == "enumerate":
        result = exact_optimize_small(
            fixture,
            state_cap_per_component=args.state_cap_per_component,
        )
        output = result.to_dict()
        output["solver"] = "exact-full-enumeration"

    elif args.solver == "layer-dp":
        dp = exact_optimize_layer_dp(
            fixture,
            orientation_cap_per_component=args.orientation_cap_per_component,
            permutation_cap_per_species=args.permutation_cap_per_species,
            transition_cap_per_component=args.transition_cap_per_component,
        )
        result = dp.layout
        output = dp.to_dict()

    elif args.solver == "branch-bound":
        bounded = optimize_branch_and_bound(
            fixture,
            node_cap_per_component=args.branch_node_cap_per_component,
            orientation_cap_per_component=args.orientation_cap_per_component,
            local_restarts=args.local_restarts,
            seed=args.seed,
        )
        result = bounded.layout
        output = bounded.to_dict()

    elif args.solver == "local-search":
        local = optimize_local_search(
            fixture,
            restarts=args.local_restarts,
            max_improving_steps=args.local_max_improving_steps,
            seed=args.seed,
        )
        result = local.layout
        output = local.to_dict()

    else:
        automatic = optimize_auto(
            fixture,
            orientation_cap_per_component=args.orientation_cap_per_component,
            permutation_cap_per_species=args.permutation_cap_per_species,
            transition_cap_per_component=args.transition_cap_per_component,
            branch_node_cap_per_component=args.branch_node_cap_per_component,
            local_restarts=args.local_restarts,
            local_max_improving_steps=args.local_max_improving_steps,
            seed=args.seed,
        )
        result = automatic.layout
        output = automatic.to_dict()

    audit = build_layout_audit(fixture, result)

    if args.audit_json:
        Path(args.audit_json).write_text(
            json.dumps(audit.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.svg:
        write_layout_comparison_svg(fixture, result, args.svg)

    output["audit"] = audit.to_dict()
    print(json.dumps(output, indent=2 if args.pretty else None, sort_keys=True))
    return 0


def entrypoint():
    import sys
    try:
        return main()
    except (ValueError, FileNotFoundError, ImportError) as error:
        print(f'syntangle: {error}',file=sys.stderr)
        if isinstance(error, ImportError):
            print('Install the optimization and plotting extras for the full workflow.',file=sys.stderr)
        return 2
