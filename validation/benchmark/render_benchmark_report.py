from __future__ import annotations

import argparse
import csv
import html
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean


MODE_ORDER = ("mild", "strong", "random")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def integer(row: dict[str, str], key: str, default: int = 0) -> int:
    value = row.get(key, "")
    if value == "":
        return default
    return int(float(value))


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def svg_shell(width: int, height: int, body: str, title: str) -> str:
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
        f'height="{height}" viewBox="0 0 {width} {height}" '
        'role="img" aria-label="' + esc(title) + '">'
        '<style>'
        'text{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;'
        'fill:#1f2937} .axis{stroke:#6b7280;stroke-width:1}'
        '.grid{stroke:#e5e7eb;stroke-width:1}.muted{fill:#6b7280}'
        '.exact{fill:#2563eb}.bounded{fill:#d97706}.initial{fill:#9ca3af}'
        '.optimum{fill:#2563eb}.mark{stroke:#111827;stroke-width:1.3;fill:#fff}'
        '</style>'
        f'{body}</svg>'
    )


def runtime_scaling_svg(rows: list[dict[str, str]]) -> str:
    width, height = 920, 520
    left, right, top, bottom = 74, 24, 48, 72
    plot_w = width - left - right
    plot_h = height - top - bottom

    x_values = [integer(row, "ancestor_chromosomes") for row in rows]
    y_values = [number(row, "wall_seconds") for row in rows]
    x_min = min(x_values)
    x_max = max(x_values)
    y_max = max(y_values) if y_values else 1.0
    y_max = max(y_max, 0.001)

    def sx(value: float) -> float:
        if x_max == x_min:
            return left + plot_w / 2
        return left + (value - x_min) / (x_max - x_min) * plot_w

    def sy(value: float) -> float:
        # log1p keeps sub-second and slower cases visible together.
        return top + plot_h - math.log1p(value) / math.log1p(y_max) * plot_h

    parts = [
        f'<text x="{left}" y="26" font-size="20" font-weight="700">'
        'Runtime scaling</text>',
        f'<text x="{left}" y="43" font-size="12" class="muted">'
        'Wall time versus ancestral chromosome count; vertical scale is log(1 + seconds).'
        '</text>',
    ]

    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        value = math.expm1(math.log1p(y_max) * frac)
        y = sy(value)
        parts.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" '
            f'y2="{y:.2f}" class="grid"/>'
        )
        parts.append(
            f'<text x="{left - 9}" y="{y + 4:.2f}" text-anchor="end" '
            f'font-size="11">{value:.2g}</text>'
        )

    for x in sorted(set(x_values)):
        px = sx(x)
        parts.append(
            f'<line x1="{px:.2f}" y1="{top}" x2="{px:.2f}" '
            f'y2="{top + plot_h}" class="grid"/>'
        )
        parts.append(
            f'<text x="{px:.2f}" y="{top + plot_h + 20}" '
            f'text-anchor="middle" font-size="11">{x}</text>'
        )

    parts.extend([
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" '
        f'y2="{top + plot_h}" class="axis"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" '
        f'y2="{top + plot_h}" class="axis"/>',
        f'<text x="{left + plot_w / 2}" y="{height - 22}" '
        f'text-anchor="middle" font-size="12">Ancestral chromosomes</text>',
        f'<text x="18" y="{top + plot_h / 2}" text-anchor="middle" '
        f'font-size="12" transform="rotate(-90 18 {top + plot_h / 2})">'
        'Wall seconds</text>',
    ])

    for row in rows:
        x = sx(integer(row, "ancestor_chromosomes"))
        y = sy(number(row, "wall_seconds"))
        is_branch = "branch-and-bound" in row.get("solver", "")
        if is_branch:
            parts.append(
                f'<rect x="{x - 4:.2f}" y="{y - 4:.2f}" width="8" height="8" '
                f'class="bounded"><title>{esc(row["case_id"])}: '
                f'{number(row, "wall_seconds"):.4f}s</title></rect>'
            )
        else:
            parts.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" class="exact">'
                f'<title>{esc(row["case_id"])}: '
                f'{number(row, "wall_seconds"):.4f}s</title></circle>'
            )

    parts.append(
        f'<text x="{left + 8}" y="{height - 46}" font-size="11">'
        'circle = exact layer DP; square = monotone component branch-and-bound'
        '</text>'
    )
    return svg_shell(width, height, "".join(parts), "Runtime scaling")


def paired_presentation_svg(rows: list[dict[str, str]]) -> str | None:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        biology_id = row.get("biology_id", "")
        if biology_id:
            grouped[biology_id].append(row)

    groups = [
        (biology_id, group)
        for biology_id, group in sorted(grouped.items())
        if len(group) >= 3
        and {item.get("tangle_mode") for item in group}
        >= set(MODE_ORDER)
    ]
    if not groups:
        return None

    width = 980
    row_h = 48
    top, left, right, bottom = 60, 190, 36, 44
    height = top + len(groups) * row_h + bottom
    plot_w = width - left - right
    max_c = max(
        integer(row, "initial_crossings")
        for _, group in groups
        for row in group
    )
    max_c = max(max_c, 1)

    def sx(value: int) -> float:
        return left + value / max_c * plot_w

    parts = [
        f'<text x="{left}" y="26" font-size="20" font-weight="700">'
        'Controlled presentation-tangle experiment</text>',
        f'<text x="{left}" y="44" font-size="12" class="muted">'
        'Each row is one biological simulation shown three different ways. '
        'Grey points are input C; blue diamonds are proven C*.'
        '</text>',
    ]

    for tick in range(5):
        value = round(max_c * tick / 4)
        x = sx(value)
        parts.append(
            f'<line x1="{x:.2f}" y1="{top - 8}" x2="{x:.2f}" '
            f'y2="{height - bottom}" class="grid"/>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{height - 16}" text-anchor="middle" '
            f'font-size="11">{value}</text>'
        )

    offsets = {"mild": -11, "strong": 0, "random": 11}
    for idx, (biology_id, group) in enumerate(groups):
        y0 = top + idx * row_h + row_h / 2
        parts.append(
            f'<text x="{left - 12}" y="{y0 + 4:.2f}" text-anchor="end" '
            f'font-size="11">{esc(biology_id)}</text>'
        )
        by_mode = {row["tangle_mode"]: row for row in group}
        optima = {
            integer(row, "optimized_crossings")
            for row in group
            if row.get("optimality_status") == "proven optimum"
        }
        for mode in MODE_ORDER:
            row = by_mode[mode]
            y = y0 + offsets[mode]
            initial = integer(row, "initial_crossings")
            x = sx(initial)
            parts.append(
                f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" class="initial">'
                f'<title>{mode}: input C={initial}</title></circle>'
            )
        if len(optima) == 1:
            optimum = next(iter(optima))
            x = sx(optimum)
            points = (
                f'{x:.2f},{y0 - 6:.2f} {x + 6:.2f},{y0:.2f} '
                f'{x:.2f},{y0 + 6:.2f} {x - 6:.2f},{y0:.2f}'
            )
            parts.append(
                f'<polygon points="{points}" class="optimum">'
                f'<title>proven C*={optimum}</title></polygon>'
            )

    return svg_shell(
        width,
        height,
        "".join(parts),
        "Controlled presentation-tangle experiment",
    )


def complexity_ladder_svg(rows: list[dict[str, str]]) -> str | None:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        biology_id = row.get("biology_id", "").strip()
        if biology_id and row.get("events_per_branch", "") != "":
            grouped[biology_id].append(row)

    groups = []
    for biology_id, group in grouped.items():
        modes = {row.get("tangle_mode") for row in group}
        if modes >= set(MODE_ORDER):
            sample = group[0]
            groups.append(
                (
                    integer(sample, "species_count"),
                    integer(sample, "ancestor_chromosomes"),
                    integer(sample, "events_per_branch"),
                    biology_id,
                    group,
                )
            )

    if len(groups) < 2:
        return None

    groups.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    width, height = 940, 500
    left, right, top, bottom = 78, 28, 54, 95
    plot_w = width - left - right
    plot_h = height - top - bottom

    ymax = max(
        number(row, "wall_seconds")
        for _, _, _, _, group in groups
        for row in group
    )
    ymax = max(ymax, 0.001)

    def sx(index: int) -> float:
        if len(groups) == 1:
            return left + plot_w / 2
        return left + index / (len(groups) - 1) * plot_w

    def sy(value: float) -> float:
        return top + plot_h - math.log1p(value) / math.log1p(ymax) * plot_h

    parts = [
        f'<text x="{left}" y="27" font-size="20" font-weight="700">'
        'Coupled complexity ladder</text>',
        f'<text x="{left}" y="45" font-size="12" class="muted">'
        'Species, chromosomes, and structural events per lineage step increase together; '
        'vertical scale is log(1 + seconds).</text>',
    ]

    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        value = math.expm1(math.log1p(ymax) * frac)
        y = sy(value)
        parts.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" '
            f'y2="{y:.2f}" class="grid"/>'
        )
        parts.append(
            f'<text x="{left - 9}" y="{y + 4:.2f}" text-anchor="end" '
            f'font-size="11">{value:.2g}</text>'
        )

    mode_classes = {"mild": "exact", "strong": "bounded", "random": "mark"}
    for mode in MODE_ORDER:
        points = []
        for index, (species, chroms, events, biology_id, group) in enumerate(groups):
            row = next(item for item in group if item.get("tangle_mode") == mode)
            x = sx(index)
            y = sy(number(row, "wall_seconds"))
            points.append((x, y, row))
        if len(points) >= 2:
            parts.append(
                '<polyline fill="none" stroke="#6b7280" stroke-width="1.5" '
                'stroke-opacity="0.55" points="{}"/>'.format(
                    " ".join(f"{x:.2f},{y:.2f}" for x, y, _ in points)
                )
            )
        for x, y, row in points:
            cls = mode_classes[mode]
            if mode == "random":
                parts.append(
                    f'<circle cx="{x:.2f}" cy="{y:.2f}" r="5" class="{cls}">'
                    f'<title>{mode}: {number(row, "wall_seconds"):.4f}s</title></circle>'
                )
            else:
                parts.append(
                    f'<circle cx="{x:.2f}" cy="{y:.2f}" r="5" class="{cls}">'
                    f'<title>{mode}: {number(row, "wall_seconds"):.4f}s</title></circle>'
                )

    for index, (species, chroms, events, biology_id, group) in enumerate(groups):
        x = sx(index)
        parts.append(
            f'<line x1="{x:.2f}" y1="{top + plot_h}" x2="{x:.2f}" '
            f'y2="{top + plot_h + 5}" class="axis"/>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{top + plot_h + 23}" text-anchor="middle" '
            f'font-size="11">{species} sp</text>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{top + plot_h + 38}" text-anchor="middle" '
            f'font-size="11">{chroms} chr</text>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{top + plot_h + 53}" text-anchor="middle" '
            f'font-size="11">{events} events/branch</text>'
        )

    parts.extend([
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" '
        f'y2="{top + plot_h}" class="axis"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" '
        f'y2="{top + plot_h}" class="axis"/>',
        f'<text x="18" y="{top + plot_h / 2}" text-anchor="middle" '
        f'font-size="12" transform="rotate(-90 18 {top + plot_h / 2})">'
        'Wall seconds</text>',
        f'<text x="{left + 12}" y="{height - 18}" font-size="11">'
        'mild = blue · strong = orange · random = outlined</text>',
    ])

    return svg_shell(
        width,
        height,
        "".join(parts),
        "Coupled complexity ladder",
    )


def orientation_reduction_svg(rows: list[dict[str, str]]) -> str | None:
    branch_rows = [
        row
        for row in rows
        if row.get("max_implicit_orientation_states", "") != ""
    ]
    if not branch_rows:
        return None

    width, height = 920, 500
    left, right, top, bottom = 86, 28, 52, 70
    plot_w = width - left - right
    plot_h = height - top - bottom

    xs = [
        max(integer(row, "max_implicit_orientation_states"), 1)
        for row in branch_rows
    ]
    ys = [
        max(
            integer(row, "total_orientation_nodes")
            + integer(row, "total_order_nodes"),
            1,
        )
        for row in branch_rows
    ]
    xmax = max(xs)
    ymax = max(ys)

    def sx(value: int) -> float:
        return left + math.log10(value) / max(math.log10(xmax), 1.0) * plot_w

    def sy(value: int) -> float:
        return (
            top
            + plot_h
            - math.log10(value) / max(math.log10(ymax), 1.0) * plot_h
        )

    parts = [
        f'<text x="{left}" y="27" font-size="20" font-weight="700">'
        'Implicit orientation space versus evaluated residual search</text>',
        f'<text x="{left}" y="45" font-size="12" class="muted">'
        'Both axes are log10. Points far below the diagonal indicate useful monotone reduction.'
        '</text>',
    ]

    for exponent in range(0, int(math.ceil(math.log10(xmax))) + 1):
        value = 10 ** exponent
        x = sx(value)
        parts.append(
            f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" '
            f'y2="{top + plot_h}" class="grid"/>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="{top + plot_h + 20}" '
            f'text-anchor="middle" font-size="11">10^{exponent}</text>'
        )

    for exponent in range(0, int(math.ceil(math.log10(ymax))) + 1):
        value = 10 ** exponent
        y = sy(value)
        parts.append(
            f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" '
            f'y2="{y:.2f}" class="grid"/>'
        )
        parts.append(
            f'<text x="{left - 10}" y="{y + 4:.2f}" '
            f'text-anchor="end" font-size="11">10^{exponent}</text>'
        )

    diagonal_max = min(xmax, ymax)
    if diagonal_max >= 1:
        parts.append(
            f'<line x1="{sx(1):.2f}" y1="{sy(1):.2f}" '
            f'x2="{sx(diagonal_max):.2f}" y2="{sy(diagonal_max):.2f}" '
            'stroke="#9ca3af" stroke-dasharray="4 4"/>'
        )

    for row, xvalue, yvalue in zip(branch_rows, xs, ys):
        x = sx(xvalue)
        y = sy(yvalue)
        parts.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="5" class="exact">'
            f'<title>{esc(row["case_id"])}: implicit={xvalue}, '
            f'evaluated={yvalue}</title></circle>'
        )
        parts.append(
            f'<text x="{x + 7:.2f}" y="{y - 7:.2f}" font-size="9">'
            f'{esc(row["case_id"])}</text>'
        )

    parts.extend([
        f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" '
        f'y2="{top + plot_h}" class="axis"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" '
        f'y2="{top + plot_h}" class="axis"/>',
        f'<text x="{left + plot_w / 2}" y="{height - 22}" '
        f'text-anchor="middle" font-size="12">Implicit legal orientation states</text>',
        f'<text x="18" y="{top + plot_h / 2}" text-anchor="middle" '
        f'font-size="12" transform="rotate(-90 18 {top + plot_h / 2})">'
        'Orientation + order nodes evaluated</text>',
    ])

    return svg_shell(
        width,
        height,
        "".join(parts),
        "Orientation-space reduction",
    )


def stage_audit(rows: list[dict[str, str]]) -> list[tuple[str, str, str]]:
    solver_counts = Counter(row.get("solver", "") for row in rows)
    normalized_helped = sum(
        integer(row, "normalization_crossings_removed") != 0
        for row in rows
    )
    multi_component = sum(
        integer(row, "incidence_component_count") > 1
        for row in rows
    )
    branch_rows = [
        row
        for row in rows
        if "branch-and-bound" in row.get("solver", "")
    ]
    forced = sum(
        integer(row, "total_orientation_groups_forced")
        for row in branch_rows
    )
    memo = sum(
        integer(row, "total_reduction_memo_hits")
        for row in branch_rows
    )

    return [
        (
            "Hypergraph incidence representation + connected components",
            "ACTIVE",
            f"Preserves n-ary homology and exactly factors independent problems; "
            f"{multi_component}/{len(rows)} cases had multiple incidence components.",
        ),
        (
            "Component canonicalization",
            "ACTIVE",
            f"Removed presentation crossings before combinatorial search in "
            f"{normalized_helped}/{len(rows)} cases.",
        ),
        (
            "GF(2) orientation basis",
            "ACTIVE",
            "Hard orientation equations and free flip groups define the legal "
            "orientation search space.",
        ),
        (
            "Exact species-layer dynamic programming",
            "ACTIVE",
            f"Primary solver for {solver_counts.get('exact-layer-dynamic-programming', 0)}"
            f"/{len(rows)} cases.",
        ),
        (
            "Monotone component branch-and-bound",
            "ACTIVE",
            f"Fallback for {len(branch_rows)}/{len(rows)} cases; independent "
            f"incidence components can run in separate worker processes. Bound "
            f"propagation recorded {forced} forcing events and {memo} memo hits.",
        ),
        (
            "Residual variable/factor graph",
            "ACTIVE DIAGNOSTIC / NEXT FACTORIZATION LAYER",
            "Records unresolved order/orientation variables, exact crossing-factor "
            "couplings, disconnected residual pieces, articulation variables, and "
            "a min-fill treewidth upper bound.",
        ),
        (
            "One-layer subset DP",
            "ACTIVE",
            "Used inside local-search incumbents and branch-and-bound conditional "
            "edge bounds/order optimization.",
        ),
        (
            "Structural bundle projection",
            "DIAGNOSTIC NOW",
            "Used to measure chromosome-level topology and benchmark complexity, "
            "but it does not currently restrict the solver search.",
        ),
        (
            "Bridges / articulation points / 2-core / biconnected kernels",
            "DIAGNOSTIC NOW",
            "Computed and visualized, but the solver does not assume these raw "
            "structural kernels factor the residual crossing objective.",
        ),
        (
            "Fundamental cycle basis",
            "DIAGNOSTIC NOW",
            "Useful structural descriptor; not currently part of the optimization "
            "objective or branching logic.",
        ),
        (
            "Derived interchromosomal ordering constraints",
            "INFRASTRUCTURE",
            "Implemented and tested, but not on the current optimizer hot path "
            "because benchmark inputs do not yet provide genuine precedence or "
            "consecutive-ones constraints.",
        ),
        (
            "Spectral ordering / braid-inspired moves",
            "DEFERRED",
            "Not used by the current solver.",
        ),
    ]


def render_html(rows: list[dict[str, str]], title: str) -> tuple[str, dict[str, str]]:
    runtime_svg = runtime_scaling_svg(rows)
    paired_svg = paired_presentation_svg(rows)
    complexity_svg = complexity_ladder_svg(rows)
    reduction_svg = orientation_reduction_svg(rows)
    svgs = {"runtime_scaling.svg": runtime_svg}
    if paired_svg is not None:
        svgs["paired_presentation.svg"] = paired_svg
    if complexity_svg is not None:
        svgs["complexity_ladder.svg"] = complexity_svg
    if reduction_svg is not None:
        svgs["orientation_reduction.svg"] = reduction_svg

    proven = sum(
        row.get("optimality_status") == "proven optimum"
        for row in rows
    )
    mean_runtime = mean(number(row, "wall_seconds") for row in rows)
    removed = sum(integer(row, "crossings_removed") for row in rows)
    initial = sum(integer(row, "initial_crossings") for row in rows)
    removed_pct = 100.0 * removed / initial if initial else 0.0

    audit_rows = "".join(
        "<tr>"
        f"<td>{esc(stage)}</td><td><strong>{esc(status)}</strong></td>"
        f"<td>{esc(note)}</td></tr>"
        for stage, status, note in stage_audit(rows)
    )

    case_rows = "".join(
        "<tr>"
        f"<td>{esc(row['case_id'])}</td>"
        f"<td>{esc(row.get('tangle_mode', ''))}</td>"
        f"<td>{integer(row, 'species_count')}</td>"
        f"<td>{integer(row, 'ancestor_chromosomes')}</td>"
        f"<td>{integer(row, 'events_per_branch')}</td>"
        f"<td>{esc(row.get('solver', ''))}</td>"
        f"<td>{integer(row, 'initial_crossings')}</td>"
        f"<td>{integer(row, 'optimized_crossings')}</td>"
        f"<td>{integer(row, 'optimality_gap')}</td>"
        f"<td>{number(row, 'wall_seconds'):.4f}</td>"
        "</tr>"
        for row in rows
    )

    figure_html = "".join(
        f'<section class="figure">{svg}</section>'
        for svg in svgs.values()
    )

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{esc(title)}</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
margin:0;background:#f8fafc;color:#111827}}
main{{max-width:1120px;margin:0 auto;padding:32px 24px 72px}}
h1{{margin-bottom:6px}} h2{{margin-top:38px}}
.sub{{color:#6b7280;margin-top:0}}
.cards{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}}
.card,.figure,table{{background:white;border:1px solid #e5e7eb;border-radius:12px}}
.card{{padding:16px}} .big{{font-size:28px;font-weight:700}}
.figure{{padding:12px;margin:16px 0;overflow-x:auto}}
svg{{max-width:100%;height:auto}}
table{{width:100%;border-collapse:separate;border-spacing:0;overflow:hidden}}
th,td{{padding:9px 10px;border-bottom:1px solid #e5e7eb;text-align:left;
font-size:13px;vertical-align:top}}
th{{background:#f3f4f6}} tr:last-child td{{border-bottom:0}}
code{{background:#eef2ff;padding:1px 4px;border-radius:4px}}
@media(max-width:800px){{.cards{{grid-template-columns:1fr 1fr}}}}
</style>
</head>
<body><main>
<h1>{esc(title)}</h1>
<p class="sub">Vector benchmark diagnostics generated directly from the TSV results.</p>
<div class="cards">
<div class="card"><div class="big">{len(rows)}</div><div>cases</div></div>
<div class="card"><div class="big">{proven}/{len(rows)}</div><div>proven optima</div></div>
<div class="card"><div class="big">{mean_runtime:.2f}s</div><div>mean wall time</div></div>
<div class="card"><div class="big">{removed_pct:.1f}%</div><div>input crossings removed</div></div>
</div>
<h2>Visual diagnostics</h2>
{figure_html}
<h2>Method-stage utility audit</h2>
<p>This separates stages that currently influence the optimizer from graph-theory
and ordering machinery that is presently diagnostic or preparatory. A stage marked
diagnostic is not being claimed as computationally useful until an ablation or
integration test demonstrates that it changes search cost or correctness.</p>
<table><thead><tr><th>Stage</th><th>Status</th><th>Current role</th></tr></thead>
<tbody>{audit_rows}</tbody></table>
<h2>Case results</h2>
<table><thead><tr><th>Case</th><th>Tangle</th><th>Species</th>
<th>Ancestor chr</th><th>Events/branch</th><th>Solver</th>
<th>Input C</th><th>Final C</th><th>Gap</th><th>Seconds</th></tr></thead>
<tbody>{case_rows}</tbody></table>
</main></body></html>
"""
    return document, svgs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_tsv")
    parser.add_argument("output_html")
    parser.add_argument("--title", default="SynTangle benchmark report")
    parser.add_argument(
        "--svg-dir",
        help="Optional directory for standalone vector SVG copies.",
    )
    args = parser.parse_args()

    rows = read_tsv(Path(args.results_tsv))
    if not rows:
        raise SystemExit("Results table is empty")

    document, svgs = render_html(rows, args.title)
    output = Path(args.output_html)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")

    if args.svg_dir:
        svg_dir = Path(args.svg_dir)
        svg_dir.mkdir(parents=True, exist_ok=True)
        for name, svg in svgs.items():
            (svg_dir / name).write_text(svg + "\n", encoding="utf-8")

    print(output)
    for name in svgs:
        print(name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
