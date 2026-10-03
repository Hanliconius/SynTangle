from __future__ import annotations

import argparse
import json
from pathlib import Path

from .audit import build_layout_audit
from .decomposition import decompose_incidence_graph
from .fixtures import load_fixture
from .incidence import analyze_fixture, build_incidence_graph
from .layer_dp import exact_optimize_layer_dp
from .layout import exact_optimize_small
from .ordering import derive_ordering_constraints
from .orientation import solve_orientation_constraints
from .visualize import write_layout_comparison_svg


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="syntangle")
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser(
        "analyze-fixture",
        help="parse a synthetic fixture and summarize its constrained structure",
    )
    analyze.add_argument("fixture")
    analyze.add_argument("--pretty", action="store_true")

    optimize = subparsers.add_parser(
        "optimize-fixture",
        help="prove a minimum crossing layout for an unambiguous fixture",
    )
    optimize.add_argument("fixture")
    optimize.add_argument(
        "--solver",
        choices=("layer-dp", "enumerate"),
        default="layer-dp",
    )
    optimize.add_argument("--state-cap-per-component", type=int, default=250000)
    optimize.add_argument("--orientation-cap-per-component", type=int, default=4096)
    optimize.add_argument("--permutation-cap-per-species", type=int, default=40320)
    optimize.add_argument("--transition-cap-per-component", type=int, default=5000000)
    optimize.add_argument("--audit-json")
    optimize.add_argument("--svg")
    optimize.add_argument("--pretty", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "analyze-fixture":
        fixture = load_fixture(args.fixture)
        graph = build_incidence_graph(fixture)
        output = analyze_fixture(fixture).to_dict()
        output["orientation"] = solve_orientation_constraints(fixture).to_dict()
        output["decomposition"] = decompose_incidence_graph(graph).to_dict()
        output["ordering_constraints"] = derive_ordering_constraints(fixture).to_dict()
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
    else:
        dp = exact_optimize_layer_dp(
            fixture,
            orientation_cap_per_component=args.orientation_cap_per_component,
            permutation_cap_per_species=args.permutation_cap_per_species,
            transition_cap_per_component=args.transition_cap_per_component,
        )
        result = dp.layout
        output = dp.to_dict()

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
