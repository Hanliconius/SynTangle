from __future__ import annotations

import argparse
import csv
import html
from collections import defaultdict
from pathlib import Path

from syntangle import (
    load_validation_bundle,
    optimize_auto,
    write_validation_report_html,
)

from generate_case_report import load_hidden_native_state


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def render_case(
    root: Path,
    row: dict[str, str],
    output_dir: Path,
    *,
    transition_cap: int,
    branch_node_cap: int,
    local_restarts: int,
) -> tuple[str, str, int, int, str]:
    case_dir = root / row["case_dir"]
    fixture = load_validation_bundle(case_dir)

    result = optimize_auto(
        fixture,
        transition_cap_per_component=transition_cap,
        branch_node_cap_per_component=branch_node_cap,
        local_restarts=local_restarts,
        seed=int(row["seed"]),
    )

    # As in generate_case_report.py, reveal hidden native state only after solve.
    native = load_hidden_native_state(fixture, case_dir)

    filename = f"{row['case_id']}_report.html"
    output = output_dir / filename
    write_validation_report_html(
        fixture,
        result.layout,
        str(output),
        solver=result.solver,
        solver_details=result.details,
        hidden_native_state=native,
    )
    return (
        row["case_id"],
        filename,
        result.layout.initial_score.crossings,
        result.layout.optimized_score.crossings,
        result.layout.optimality_status,
    )


def write_index(
    root: Path,
    rendered: list[tuple[dict[str, str], tuple[str, str, int, int, str]]],
    output_dir: Path,
) -> None:
    grouped = defaultdict(list)
    for row, result in rendered:
        grouped[row.get("biology_id", row["case_id"])].append((row, result))

    sections = []
    for biology_id, items in sorted(grouped.items()):
        sample = items[0][0]
        meta = (
            f"{sample.get('species_count', '?')} species · "
            f"{sample.get('ancestor_chromosomes', '?')} ancestral chromosomes · "
            f"{sample.get('events_per_branch', '?')} events/branch · "
            f"{sample.get('genes_per_chromosome', '?')} anchors/chromosome"
        )
        cards = []
        for row, result in sorted(
            items,
            key=lambda item: ("mild", "strong", "random").index(
                item[0].get("tangle_mode", "random")
            ),
        ):
            case_id, filename, initial, final, status = result
            cards.append(
                f"""
                <a class="card" href="{html.escape(filename)}">
                  <div class="mode">{html.escape(row.get('tangle_mode', ''))}</div>
                  <div class="score">{initial} → {final}</div>
                  <div class="status">{html.escape(status)}</div>
                  <div class="open">Open full native / tangled / optimized audit →</div>
                </a>
                """
            )
        sections.append(
            f"""
            <section>
              <h2>{html.escape(biology_id)}</h2>
              <p>{html.escape(meta)}</p>
              <div class="cards">{''.join(cards)}</div>
            </section>
            """
        )

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SynTangle visual ground-truth gallery</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
background:#f8fafc;color:#111827;margin:0}}
main{{max-width:1100px;margin:0 auto;padding:32px 24px 72px}}
h1{{margin-bottom:6px}} .sub{{color:#6b7280}}
section{{margin-top:34px}} section>p{{color:#4b5563}}
.cards{{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}}
.card{{display:block;text-decoration:none;color:inherit;background:white;
border:1px solid #e5e7eb;border-radius:12px;padding:18px}}
.card:hover{{border-color:#9ca3af}}
.mode{{font-size:13px;text-transform:uppercase;letter-spacing:.05em;color:#6b7280}}
.score{{font-size:28px;font-weight:700;margin:8px 0}}
.status{{font-size:13px;color:#374151}}
.open{{font-size:13px;margin-top:14px;color:#1d4ed8}}
@media(max-width:760px){{.cards{{grid-template-columns:1fr}}}}
</style>
</head>
<body><main>
<h1>SynTangle visual ground-truth gallery</h1>
<p class="sub">
Each card opens the full case audit: hidden simulator-native baseline,
deliberately tangled public input, optimized output, graph structure, and
solver details.
</p>
{''.join(sections)}
</main></body></html>
"""
    (output_dir / "index.html").write_text(document, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark_root")
    parser.add_argument("output_dir")
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=("mild", "strong", "random"),
        default=("mild", "strong", "random"),
    )
    parser.add_argument("--transition-cap", type=int, default=250000)
    parser.add_argument("--branch-node-cap", type=int, default=100000)
    parser.add_argument("--local-restarts", type=int, default=4)
    args = parser.parse_args()

    root = Path(args.benchmark_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = read_tsv(root / "benchmark_manifest.tsv")
    rows = [
        row for row in manifest
        if row.get("tangle_mode") in set(args.modes)
    ]
    if not rows:
        raise SystemExit("No benchmark cases matched the requested modes")

    rendered = []
    for row in rows:
        result = render_case(
            root,
            row,
            output_dir,
            transition_cap=args.transition_cap,
            branch_node_cap=args.branch_node_cap,
            local_restarts=args.local_restarts,
        )
        rendered.append((row, result))
        print(result[0], result[2], "->", result[3], result[4])

    write_index(root, rendered, output_dir)
    print(output_dir / "index.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
