from __future__ import annotations

import argparse
import json

from .decomposition import decompose_incidence_graph
from .fixtures import load_fixture
from .incidence import analyze_fixture, build_incidence_graph
from .layout import exact_optimize_small
from .ordering import derive_ordering_constraints
from .orientation import solve_orientation_constraints


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
        help="prove a minimum crossing layout for a small unambiguous fixture",
    )
    optimize.add_argument("fixture")
    optimize.add_argument("--state-cap-per-component", type=int, default=250000)
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

    if args.command == "optimize-fixture":
        fixture = load_fixture(args.fixture)
        result = exact_optimize_small(
            fixture,
            state_cap_per_component=args.state_cap_per_component,
        )
        print(json.dumps(result.to_dict(), indent=2 if args.pretty else None, sort_keys=True))
        return 0

    raise AssertionError(f"Unhandled command: {args.command}")
