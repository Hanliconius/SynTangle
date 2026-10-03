from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from statistics import mean
from time import perf_counter

from syntangle import (
    LayoutState,
    build_incidence_graph,
    build_residual_factorization,
    build_structural_projection,
    fixture_fingerprint,
    load_validation_bundle,
    optimize_auto,
    score_crossings,
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def raw_cycle_rank(fixture) -> int:
    graph = build_incidence_graph(fixture)
    return sum(
        graph.summarize_component(component).cycle_rank
        for component in graph.connected_components()
    )


def biological_fingerprint(fixture) -> str:
    """Hash extant biological evidence while ignoring presentation state."""

    parts: list[str] = []
    for chromosome in sorted(fixture.chromosomes, key=lambda item: item.ref):
        parts.append(
            "C\t{}\t{}\t{:.12g}".format(
                chromosome.ref.species_id,
                chromosome.ref.chromosome_id,
                chromosome.length,
            )
        )
        for block in chromosome.blocks:
            parts.append(
                "B\t{}\t{}\t{}\t{:.12g}\t{:.12g}\t{}".format(
                    chromosome.ref.species_id,
                    chromosome.ref.chromosome_id,
                    block.homology_id,
                    block.start,
                    block.end,
                    block.strand,
                )
            )
    for constraint in sorted(
        fixture.orientation_constraints,
        key=lambda item: (item.a, item.b, item.xor),
    ):
        parts.append(
            "O\t{}\t{}\t{}".format(
                constraint.a.label,
                constraint.b.label,
                constraint.xor,
            )
        )
    payload = "\n".join(parts).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def hidden_native_score(fixture, case_dir: Path) -> int | None:
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
            orientation[ref_lookup[(species, row["chrom"])]] = int(
                row["orientation"]
            )

    state = LayoutState(
        chromosome_order=order,
        chromosome_orientation=orientation,
    )
    return score_crossings(fixture, state).crossings


def hidden_counts(case_dir: Path) -> tuple[int, int, int]:
    evolution_path = case_dir / "hidden_evolution_log.tsv"
    tangle_path = case_dir / "hidden_tangle_log.tsv"

    event_count = 0
    event_types: set[str] = set()
    if evolution_path.is_file():
        rows = read_tsv(evolution_path)
        event_count = len(rows)
        event_types = {
            row.get("type", "")
            for row in rows
            if row.get("type")
        }

    changed_layout_rows = 0
    if tangle_path.is_file():
        rows = read_tsv(tangle_path)
        for row in rows:
            rank_changed = (
                row.get("display_rank_before")
                != row.get("display_rank_after")
            )
            flipped = row.get("flipped", "").strip().upper() == "TRUE"
            if rank_changed or flipped:
                changed_layout_rows += 1

    return event_count, len(event_types), changed_layout_rows


def bounds_from_result(result) -> tuple[int, int, int]:
    optimized = result.layout.optimized_score.crossings
    if result.layout.optimality_status == "proven optimum":
        return optimized, optimized, 0

    lower = int(result.details.get("lower_bound", 0))
    upper = int(result.details.get("upper_bound", optimized))
    return lower, upper, upper - lower


def reduction_metrics(result) -> dict[str, object]:
    diagnostics = result.details.get("component_diagnostics", [])
    monotone = [
        item
        for item in diagnostics
        if "implicit_orientation_states" in item
    ]
    if not monotone:
        return {
            "max_implicit_orientation_states": "",
            "total_orientation_nodes": "",
            "total_order_nodes": "",
            "total_orientation_branches_pruned": "",
            "total_orientation_groups_forced": "",
            "total_reduction_memo_hits": "",
        }

    return {
        "max_implicit_orientation_states": max(
            int(item["implicit_orientation_states"])
            for item in monotone
        ),
        "total_orientation_nodes": sum(
            int(item["orientation_nodes_evaluated"])
            for item in monotone
        ),
        "total_order_nodes": sum(
            int(item["order_nodes_evaluated"])
            for item in monotone
        ),
        "total_orientation_branches_pruned": sum(
            int(item["orientation_branches_pruned"])
            for item in monotone
        ),
        "total_orientation_groups_forced": sum(
            int(item["orientation_groups_forced"])
            for item in monotone
        ),
        "total_reduction_memo_hits": sum(
            int(item["memo_hits"])
            for item in monotone
        ),
    }


def run_case(
    row: dict[str, str],
    root: Path,
    *,
    transition_cap: int,
    branch_node_cap: int,
    local_restarts: int,
    component_workers: int,
) -> dict[str, object]:
    case_dir = root / row["case_dir"]
    fixture = load_validation_bundle(case_dir)

    graph = build_incidence_graph(fixture)
    components = graph.connected_components()
    structural = build_structural_projection(fixture)
    residual = build_residual_factorization(fixture)
    raw_rank = sum(
        graph.summarize_component(component).cycle_rank
        for component in components
    )
    max_kernel_nodes = max(
        (len(kernel) for kernel in structural.decomposition.hard_kernels),
        default=0,
    )

    started = perf_counter()
    result = optimize_auto(
        fixture,
        transition_cap_per_component=transition_cap,
        branch_node_cap_per_component=branch_node_cap,
        local_restarts=local_restarts,
        seed=int(row["seed"]),
        component_workers=component_workers,
    )
    elapsed = perf_counter() - started

    layout = result.layout
    lower, upper, gap = bounds_from_result(result)
    reduction = reduction_metrics(result)

    if layout.optimized_score.crossings > layout.initial_score.crossings:
        raise AssertionError(
            f"{fixture.fixture_id}: optimized layout is worse than input"
        )
    if not (lower <= layout.optimized_score.crossings <= upper):
        raise AssertionError(
            f"{fixture.fixture_id}: invalid lower/upper bounds"
        )
    if layout.optimality_status == "proven optimum":
        if not (
            lower
            == upper
            == layout.optimized_score.crossings
        ):
            raise AssertionError(
                f"{fixture.fixture_id}: proven result has nonzero gap"
            )

    # Hidden simulator truth is intentionally read only after optimization.
    native_crossings = hidden_native_score(fixture, case_dir)
    event_count, event_type_count, changed_layout_rows = hidden_counts(
        case_dir
    )

    return {
        **row,
        "input_fingerprint": fixture_fingerprint(fixture),
        "biological_fingerprint": biological_fingerprint(fixture),
        "incidence_component_count": len(components),
        "chromosome_count": len(fixture.chromosomes),
        "homology_group_count": len(fixture.homology_ids),
        "structural_bundle_count": len(structural.bundles),
        "raw_cycle_rank": raw_rank,
        "structural_cycle_rank": structural.cycle_rank_total,
        "structural_hard_kernel_count": len(
            structural.decomposition.hard_kernels
        ),
        "max_structural_kernel_nodes": max_kernel_nodes,
        "residual_variable_count": len(residual.variables),
        "residual_factor_count": len(residual.factors),
        "residual_objective_component_count": (
            residual.objective_component_count
        ),
        "max_residual_component_variables": (
            residual.max_objective_component_variables
        ),
        "residual_isolated_variable_count": len(
            residual.isolated_variable_ids
        ),
        "residual_articulation_variable_count": len(
            residual.articulation_variable_ids
        ),
        "residual_treewidth_upper_bound": (
            residual.min_fill_treewidth_upper_bound
        ),
        "hidden_evolution_event_count": event_count,
        "hidden_evolution_event_type_count": event_type_count,
        "hidden_tangle_changed_chromosomes": changed_layout_rows,
        "hidden_native_crossings": (
            "" if native_crossings is None else native_crossings
        ),
        "initial_crossings": layout.initial_score.crossings,
        "normalized_crossings": layout.normalized_score.crossings,
        "normalization_crossings_removed": (
            layout.initial_score.crossings
            - layout.normalized_score.crossings
        ),
        "optimized_crossings": layout.optimized_score.crossings,
        "crossings_removed": (
            layout.initial_score.crossings
            - layout.optimized_score.crossings
        ),
        "tangle_added_crossings_vs_native": (
            ""
            if native_crossings is None
            else layout.initial_score.crossings - native_crossings
        ),
        "optimized_minus_native_crossings": (
            ""
            if native_crossings is None
            else layout.optimized_score.crossings - native_crossings
        ),
        "lower_bound": lower,
        "upper_bound": upper,
        "optimality_gap": gap,
        "optimality_status": layout.optimality_status,
        "solver": result.solver,
        "component_workers_used": int(
            result.details.get("component_workers", 1)
        ),
        "states_or_nodes_evaluated": layout.states_evaluated,
        **reduction,
        "wall_seconds": round(elapsed, 6),
    }


def write_results(rows: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    rows: list[dict[str, object]],
    json_path: Path,
    markdown_path: Path,
) -> None:
    solvers = Counter(str(row["solver"]) for row in rows)
    proven = sum(
        row["optimality_status"] == "proven optimum"
        for row in rows
    )
    summary = {
        "case_count": len(rows),
        "proven_optimum_count": proven,
        "bounded_count": len(rows) - proven,
        "mean_initial_crossings": mean(
            float(row["initial_crossings"]) for row in rows
        ),
        "mean_optimized_crossings": mean(
            float(row["optimized_crossings"]) for row in rows
        ),
        "mean_crossings_removed": mean(
            float(row["crossings_removed"]) for row in rows
        ),
        "mean_wall_seconds": mean(
            float(row["wall_seconds"]) for row in rows
        ),
        "solver_counts": dict(sorted(solvers.items())),
    }
    json_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# SynTangle simulation benchmark",
        "",
        f"- Cases: {summary['case_count']}",
        f"- Proven optima: {summary['proven_optimum_count']}",
        f"- Bounded/best-known: {summary['bounded_count']}",
        f"- Mean input crossings: {summary['mean_initial_crossings']:.3f}",
        f"- Mean optimized crossings: {summary['mean_optimized_crossings']:.3f}",
        f"- Mean crossings removed: {summary['mean_crossings_removed']:.3f}",
        f"- Mean wall time: {summary['mean_wall_seconds']:.4f} s",
        "",
        "| Case | Solver | Initial C | Final C | Gap | Structural μ | Max kernel | Seconds |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| {case_id} | {solver} | {initial_crossings} | "
            "{optimized_crossings} | {optimality_gap} | "
            "{structural_cycle_rank} | {max_structural_kernel_nodes} | "
            "{wall_seconds} |".format(**row)
        )
    markdown_path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark_root")
    parser.add_argument("--output", default="benchmark_results.tsv")
    parser.add_argument("--summary-json", default="benchmark_summary.json")
    parser.add_argument("--summary-md", default="benchmark_summary.md")
    parser.add_argument("--transition-cap", type=int, default=250000)
    parser.add_argument("--branch-node-cap", type=int, default=100000)
    parser.add_argument("--local-restarts", type=int, default=4)
    parser.add_argument(
        "--component-workers",
        type=int,
        default=1,
        help=(
            "Worker processes for independent incidence components inside "
            "branch-and-bound. Use 1 for the historical serial baseline."
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from completed cases already present in --output.",
    )
    args = parser.parse_args()

    if args.component_workers < 1:
        raise SystemExit("--component-workers must be at least 1")

    root = Path(args.benchmark_root)
    manifest = read_tsv(root / "benchmark_manifest.tsv")
    if not manifest:
        raise SystemExit("Benchmark manifest is empty")

    output = Path(args.output)
    summary_json = Path(args.summary_json)
    summary_md = Path(args.summary_md)

    completed: dict[str, dict[str, object]] = {}
    if args.resume and output.is_file():
        for saved in read_tsv(output):
            case_id = saved.get("case_id", "")
            if case_id:
                completed[case_id] = dict(saved)

    results: list[dict[str, object]] = []
    total = len(manifest)

    for index, row in enumerate(manifest, start=1):
        case_id = row["case_id"]
        if case_id in completed:
            result = completed[case_id]
            results.append(result)
            print(
                f"[{index}/{total}] SKIP {case_id} "
                "(checkpoint already complete)",
                flush=True,
            )
            continue

        print(
            f"[{index}/{total}] START {case_id}",
            flush=True,
        )
        result = run_case(
            row,
            root,
            transition_cap=args.transition_cap,
            branch_node_cap=args.branch_node_cap,
            local_restarts=args.local_restarts,
            component_workers=args.component_workers,
        )
        results.append(result)

        # Checkpoint immediately. If the next case is interrupted, every
        # completed case remains available for --resume.
        write_results(results, output)
        write_summary(results, summary_json, summary_md)

        print(
            "[{}/{}] DONE  {}  C: {} -> {}  gap={}  {}  {:.3f}s".format(
                index,
                total,
                case_id,
                result["initial_crossings"],
                result["optimized_crossings"],
                result["optimality_gap"],
                result["optimality_status"],
                float(result["wall_seconds"]),
            ),
            flush=True,
        )

    # Re-write once in manifest order after a resumed run.
    write_results(results, output)
    write_summary(results, summary_json, summary_md)

    print(summary_md.read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
