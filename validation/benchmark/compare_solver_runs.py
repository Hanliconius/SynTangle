from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["case_id"]: row for row in rows}


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in ("", None) else 0.0


def compare_rows(
    baseline: dict[str, dict[str, str]],
    candidate: dict[str, dict[str, str]],
) -> tuple[list[dict[str, object]], list[str]]:
    common = sorted(set(baseline) & set(candidate))
    output: list[dict[str, object]] = []
    errors: list[str] = []

    for case_id in common:
        before = baseline[case_id]
        after = candidate[case_id]
        before_wall = as_float(before, "wall_seconds")
        after_wall = as_float(after, "wall_seconds")
        speedup = (
            before_wall / after_wall
            if after_wall > 0
            else float("inf")
        )

        before_c = int(float(before["optimized_crossings"]))
        after_c = int(float(after["optimized_crossings"]))
        both_proven = (
            before.get("optimality_status") == "proven optimum"
            and after.get("optimality_status") == "proven optimum"
        )
        if both_proven and before_c != after_c:
            errors.append(
                f"{case_id}: proven optimum changed "
                f"{before_c} -> {after_c}"
            )

        output.append(
            {
                "case_id": case_id,
                "baseline_solver": before.get("solver", ""),
                "candidate_solver": after.get("solver", ""),
                "baseline_wall": before_wall,
                "candidate_wall": after_wall,
                "speedup": speedup,
                "baseline_c": before_c,
                "candidate_c": after_c,
                "baseline_status": before.get(
                    "optimality_status", ""
                ),
                "candidate_status": after.get(
                    "optimality_status", ""
                ),
                "residual_table_entries": after.get(
                    "total_residual_table_entries", ""
                ),
                "leaf_eliminations": after.get(
                    "total_residual_leaf_eliminations", ""
                ),
                "articulation_conditionings": after.get(
                    "total_residual_articulation_conditionings", ""
                ),
                "dynamic_splits": after.get(
                    "total_residual_dynamic_factor_splits", ""
                ),
                "min_fill_eliminations": after.get(
                    "total_residual_min_fill_eliminations", ""
                ),
                "scope_variables_removed": after.get(
                    "total_factor_scope_variables_removed", ""
                ),
                "max_intermediate_scope": after.get(
                    "max_residual_intermediate_scope", ""
                ),
            }
        )

    return output, errors


def render_markdown(
    comparisons: list[dict[str, object]],
    errors: list[str],
) -> str:
    lines = [
        "# SynTangle solver A/B comparison",
        "",
        "| Case | Old solver | New solver | Old s | New s | Speedup | "
        "C old/new | New residual work |",
        "|---|---|---|---:|---:|---:|---:|---|",
    ]

    for row in comparisons:
        work_parts = []
        for label, key in (
            ("tables", "residual_table_entries"),
            ("leaf", "leaf_eliminations"),
            ("art", "articulation_conditionings"),
            ("split", "dynamic_splits"),
            ("minfill", "min_fill_eliminations"),
            ("scope-", "scope_variables_removed"),
            ("maxscope", "max_intermediate_scope"),
        ):
            value = row[key]
            if value not in ("", None):
                work_parts.append(f"{label}={value}")
        work = ", ".join(work_parts) if work_parts else "-"

        lines.append(
            "| {case_id} | {baseline_solver} | {candidate_solver} | "
            "{baseline_wall:.3f} | {candidate_wall:.3f} | "
            "{speedup:.2f}x | {baseline_c}/{candidate_c} | {work} |".format(
                work=work,
                **row,
            )
        )

    if comparisons:
        finite = [
            float(row["speedup"])
            for row in comparisons
            if float(row["candidate_wall"]) > 0
        ]
        if finite:
            geometric = 1.0
            for value in finite:
                geometric *= value
            geometric **= 1.0 / len(finite)
            lines.extend(
                [
                    "",
                    f"- Geometric-mean speedup: {geometric:.2f}x",
                ]
            )

    lines.extend(
        [
            f"- Proven-optimum mismatches: {len(errors)}",
            "",
        ]
    )
    if errors:
        lines.append("## Correctness failures")
        lines.append("")
        lines.extend(f"- {error}" for error in errors)
        lines.append("")

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline")
    parser.add_argument("candidate")
    parser.add_argument("--markdown", required=True)
    args = parser.parse_args()

    baseline = read_rows(Path(args.baseline))
    candidate = read_rows(Path(args.candidate))
    comparisons, errors = compare_rows(baseline, candidate)

    if not comparisons:
        raise SystemExit("No common case_id values between runs")

    markdown = render_markdown(comparisons, errors)
    output = Path(args.markdown)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(markdown, encoding="utf-8")
    print(markdown)

    if errors:
        raise SystemExit(1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
