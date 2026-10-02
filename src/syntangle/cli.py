from __future__ import annotations

import argparse
import json

from .fixtures import load_fixture
from .incidence import analyze_fixture


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="syntangle")
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser(
        "analyze-fixture",
        help="parse a synthetic fixture and summarize its chromosome↔homology incidence graph",
    )
    analyze.add_argument("fixture")
    analyze.add_argument("--pretty", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "analyze-fixture":
        fixture = load_fixture(args.fixture)
        analysis = analyze_fixture(fixture)
        print(json.dumps(analysis.to_dict(), indent=2 if args.pretty else None, sort_keys=True))
        return 0
    raise AssertionError(f"Unhandled command: {args.command}")
