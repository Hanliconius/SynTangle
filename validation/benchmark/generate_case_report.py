from __future__ import annotations

import argparse
import csv
from pathlib import Path

from syntangle import (
    LayoutState,
    load_validation_bundle,
    optimize_auto,
    write_validation_report_html,
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_hidden_native_state(fixture, case_dir: Path) -> LayoutState | None:
    path = case_dir / "hidden_native_display_state.tsv"
    if not path.is_file():
        return None

    rows = read_tsv(path)
    ref_lookup = {
        (ref.species_id, ref.chromosome_id): ref
        for ref in fixture.chromosome_refs
    }
    order = {}
    orientation = {}

    for species in fixture.species_ids:
        species_rows = [
            row for row in rows if row["species"] == species
        ]
        species_rows.sort(key=lambda row: int(row["display_rank"]))
        order[species] = tuple(
            ref_lookup[(species, row["chrom"])]
            for row in species_rows
        )
        for row in species_rows:
            ref = ref_lookup[(species, row["chrom"])]
            orientation[ref] = int(row["orientation"])

    return LayoutState(
        chromosome_order=order,
        chromosome_orientation=orientation,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir")
    parser.add_argument("output_html")
    parser.add_argument("--transition-cap", type=int, default=250000)
    parser.add_argument("--branch-node-cap", type=int, default=100000)
    parser.add_argument("--local-restarts", type=int, default=4)
    parser.add_argument("--component-workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    case_dir = Path(args.case_dir)
    fixture = load_validation_bundle(case_dir)

    # Solve strictly from public bundle content.
    result = optimize_auto(
        fixture,
        transition_cap_per_component=args.transition_cap,
        branch_node_cap_per_component=args.branch_node_cap,
        local_restarts=args.local_restarts,
        seed=args.seed,
        component_workers=args.component_workers,
    )

    # Hidden native presentation is loaded only after optimization.
    native = load_hidden_native_state(fixture, case_dir)

    output = Path(args.output_html)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_validation_report_html(
        fixture,
        result.layout,
        str(output),
        solver=result.solver,
        solver_details=result.details,
        hidden_native_state=native,
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
