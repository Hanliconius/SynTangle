from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_int(row: dict[str, str], key: str) -> int:
    value = row.get(key, "")
    return int(float(value)) if value not in ("", None) else 0


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in ("", None) else 0.0


def axis_rows(
    rows: list[dict[str, str]],
    axis: str,
) -> list[dict[str, str]]:
    independent = [
        row for row in rows
        if row.get("lineage_model") == "independent"
    ]

    if axis == "species":
        selected = [
            row for row in independent
            if as_int(row, "ancestor_chromosomes") == 16
            and as_int(row, "events_per_branch") == 1
        ]
        return sorted(selected, key=lambda row: as_int(row, "species_count"))

    if axis == "chromosomes":
        selected = [
            row for row in independent
            if as_int(row, "species_count") == 8
            and as_int(row, "events_per_branch") == 1
        ]
        return sorted(
            selected,
            key=lambda row: as_int(row, "ancestor_chromosomes"),
        )

    if axis == "rearrangements":
        selected = [
            row for row in independent
            if as_int(row, "species_count") == 8
            and as_int(row, "ancestor_chromosomes") == 16
        ]
        return sorted(
            selected,
            key=lambda row: as_int(row, "events_per_branch"),
        )

    if axis == "lepidoptera_like":
        selected = [
            row for row in independent
            if as_int(row, "ancestor_chromosomes") == 31
            and as_int(row, "events_per_branch") <= 1
        ]
        return sorted(
            selected,
            key=lambda row: (
                as_int(row, "species_count"),
                as_int(row, "events_per_branch"),
            ),
        )

    raise ValueError(f"Unknown axis: {axis}")


def level_label(row: dict[str, str], axis: str) -> str:
    if axis == "species":
        return str(as_int(row, "species_count"))
    if axis == "chromosomes":
        return str(as_int(row, "ancestor_chromosomes"))
    if axis == "rearrangements":
        return str(as_int(row, "events_per_branch"))
    return (
        f"{as_int(row, 'species_count')} sp / "
        f"{as_int(row, 'events_per_branch')} ev"
    )


def section(
    rows: list[dict[str, str]],
    axis: str,
    title: str,
) -> list[str]:
    chosen = axis_rows(rows, axis)
    if not chosen:
        return [f"## {title}", "", "_No matching completed cases._", ""]

    baseline = max(as_float(chosen[0], "wall_seconds"), 1e-12)
    lines = [
        f"## {title}",
        "",
        "| Level | Wall s | Relative wall | Proof | Solver | "
        "Residual vars (max comp.) | Treewidth UB | Nodes/states | Final C |",
        "|---:|---:|---:|---|---|---:|---:|---:|---:|",
    ]

    for row in chosen:
        wall = as_float(row, "wall_seconds")
        lines.append(
            "| {} | {:.4f} | {:.2f}x | {} | {} | {} | {} | {} | {} |".format(
                level_label(row, axis),
                wall,
                wall / baseline,
                row.get("optimality_status", ""),
                row.get("solver", ""),
                as_int(row, "max_residual_component_variables"),
                as_int(row, "residual_treewidth_upper_bound"),
                as_int(row, "states_or_nodes_evaluated"),
                as_int(row, "optimized_crossings"),
            )
        )

    max_wall = max(as_float(row, "wall_seconds") for row in chosen)
    lines.extend(
        [
            "",
            (
                f"Runtime span across this controlled axis: "
                f"{max_wall / baseline:.2f}x relative to its first level."
            ),
            "",
        ]
    )
    return lines


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("results")
    parser.add_argument(
        "--markdown",
        default="factorial_analysis.md",
    )
    args = parser.parse_args()

    rows = read_rows(Path(args.results))
    if not rows:
        raise SystemExit("Results table is empty")

    lines = [
        "# SynTangle orthogonal biological benchmark",
        "",
        "This report varies one major biological axis at a time using "
        "independent descendants from a common hidden ancestor.",
        "",
        "The runtime ratios are descriptive rather than formal scaling laws. "
        "Residual component size/treewidth diagnostics are included because "
        "they may explain runtime better than raw genome dimensions.",
        "",
    ]
    lines.extend(section(rows, "species", "Species-count axis"))
    lines.extend(section(rows, "chromosomes", "Chromosome-count axis"))
    lines.extend(
        section(rows, "rearrangements", "Rearrangement-burden axis")
    )
    lines.extend(
        section(rows, "lepidoptera_like", "Lepidoptera-like 31-chromosome checks")
    )

    completed = len(rows)
    proven = sum(
        row.get("optimality_status") == "proven optimum"
        for row in rows
    )
    bounded = completed - proven
    lines.extend(
        [
            "## Completion",
            "",
            f"- Completed cases: {completed}",
            f"- Proven optima: {proven}",
            f"- Bounded/best-known: {bounded}",
            "",
        ]
    )

    output = Path(args.markdown)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8")
    print(output.read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
