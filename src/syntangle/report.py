from __future__ import annotations

from html import escape
import json

from .audit import build_layout_audit
from .decomposition import decompose_incidence_graph
from .incidence import IncidenceGraph, build_incidence_graph
from .layout import ExactLayoutResult, LayoutState, score_crossings
from .model import Fixture
from .residual import ResidualFactorization, build_residual_factorization
from .structural import build_structural_projection
from .visualize import render_layout_comparison_svg, render_layout_state_svg


def _short_node_label(node_id: str) -> str:
    if node_id.startswith("chrom::"):
        parts = node_id.split("::", 2)
        return f"{parts[1]}:{parts[2]}"
    if node_id.startswith("homology::"):
        return node_id.split("::", 1)[1]
    return node_id


def render_incidence_graph_svg(
    graph: IncidenceGraph,
    *,
    title: str,
    core_nodes: frozenset[str] = frozenset(),
    width: int = 1000,
    max_render_nodes: int = 180,
) -> str:
    """Render a deterministic chromosome↔homology graph for visual audit."""

    chromosome_nodes = tuple(
        node for node in graph.node_ids if node.startswith("chrom::")
    )
    homology_nodes = tuple(
        node for node in graph.node_ids if node.startswith("homology::")
    )

    if graph.vertex_count > max_render_nodes:
        height = 180
        return "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
                f'viewBox="0 0 {width} {height}">',
                '<rect width="100%" height="100%" fill="white"/>',
                f'<text x="{width / 2}" y="35" text-anchor="middle" '
                f'font-family="sans-serif" font-size="18">{escape(title)}</text>',
                f'<text x="{width / 2}" y="85" text-anchor="middle" '
                f'font-family="sans-serif" font-size="14">'
                f'{graph.vertex_count} vertices, {graph.edge_count} incidence edges</text>',
                f'<text x="{width / 2}" y="115" text-anchor="middle" '
                f'font-family="sans-serif" font-size="13">'
                f'Node rendering suppressed above {max_render_nodes} vertices.</text>',
                "</svg>",
            ]
        )

    row_gap = 22
    top = 60
    height = max(
        260,
        top + row_gap * max(len(chromosome_nodes), len(homology_nodes), 1) + 40,
    )
    left_x = 260
    right_x = width - 260

    chromosome_y = {
        node: top + index * row_gap
        for index, node in enumerate(chromosome_nodes)
    }
    homology_y = {
        node: top + index * row_gap
        for index, node in enumerate(homology_nodes)
    }

    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="26" text-anchor="middle" '
        f'font-family="sans-serif" font-size="18">{escape(title)}</text>',
        f'<text x="{left_x}" y="45" text-anchor="middle" '
        f'font-family="sans-serif" font-size="12">chromosomes</text>',
        f'<text x="{right_x}" y="45" text-anchor="middle" '
        f'font-family="sans-serif" font-size="12">homology / structural bundles</text>',
    ]

    for a, b in graph.edge_endpoints:
        if a.startswith("chrom::"):
            cnode, hnode = a, b
        else:
            cnode, hnode = b, a
        body.append(
            f'<line x1="{left_x}" y1="{chromosome_y[cnode]}" '
            f'x2="{right_x}" y2="{homology_y[hnode]}" '
            f'stroke="#888" stroke-opacity="0.22" stroke-width="1"/>'
        )

    for node in chromosome_nodes:
        is_core = node in core_nodes
        stroke = "#b2182b" if is_core else "#222"
        stroke_width = 2.5 if is_core else 1.2
        y = chromosome_y[node]
        body.append(
            f'<circle cx="{left_x}" cy="{y}" r="5" fill="white" '
            f'stroke="{stroke}" stroke-width="{stroke_width}"/>'
        )
        body.append(
            f'<text x="{left_x - 10}" y="{y + 4}" text-anchor="end" '
            f'font-family="sans-serif" font-size="10">'
            f'{escape(_short_node_label(node))}</text>'
        )

    for node in homology_nodes:
        is_core = node in core_nodes
        stroke = "#b2182b" if is_core else "#222"
        stroke_width = 2.5 if is_core else 1.2
        y = homology_y[node]
        body.append(
            f'<circle cx="{right_x}" cy="{y}" r="4" fill="white" '
            f'stroke="{stroke}" stroke-width="{stroke_width}"/>'
        )
        body.append(
            f'<text x="{right_x + 10}" y="{y + 4}" text-anchor="start" '
            f'font-family="sans-serif" font-size="9">'
            f'{escape(_short_node_label(node))}</text>'
        )

    body.append("</svg>")
    return "\n".join(body)


def render_residual_factor_graph_svg(
    residual: ResidualFactorization,
    *,
    title: str = "Residual decision / factor graph",
    width: int = 1000,
    max_render_nodes: int = 180,
) -> str:
    variable_ids = tuple(
        sorted(variable.variable_id for variable in residual.variables)
    )
    factor_ids = tuple(
        sorted(
            factor.factor_id
            for factor in residual.factors
            if factor.variable_ids
        )
    )
    total_nodes = len(variable_ids) + len(factor_ids)

    if total_nodes > max_render_nodes:
        height = 190
        return "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
                f'viewBox="0 0 {width} {height}">',
                '<rect width="100%" height="100%" fill="white"/>',
                f'<text x="{width / 2}" y="35" text-anchor="middle" '
                f'font-family="sans-serif" font-size="18">{escape(title)}</text>',
                f'<text x="{width / 2}" y="80" text-anchor="middle" '
                f'font-family="sans-serif" font-size="14">'
                f'{len(variable_ids)} residual variables, {len(factor_ids)} active factors</text>',
                f'<text x="{width / 2}" y="108" text-anchor="middle" '
                f'font-family="sans-serif" font-size="13">'
                f'{residual.objective_component_count} objective components; '
                f'max {residual.max_objective_component_variables} variables</text>',
                f'<text x="{width / 2}" y="136" text-anchor="middle" '
                f'font-family="sans-serif" font-size="13">'
                f'min-fill treewidth upper bound: '
                f'{residual.min_fill_treewidth_upper_bound}</text>',
                f'<text x="{width / 2}" y="164" text-anchor="middle" '
                f'font-family="sans-serif" font-size="12">'
                f'Node rendering suppressed above {max_render_nodes} nodes.</text>',
                "</svg>",
            ]
        )

    row_gap = 22
    top = 62
    height = max(
        280,
        top + row_gap * max(len(variable_ids), len(factor_ids), 1) + 44,
    )
    left_x = 310
    right_x = width - 310

    variable_y = {
        node: top + index * row_gap
        for index, node in enumerate(variable_ids)
    }
    factor_y = {
        node: top + index * row_gap
        for index, node in enumerate(factor_ids)
    }
    factor_lookup = {
        factor.factor_id: factor
        for factor in residual.factors
        if factor.variable_ids
    }
    articulation = set(residual.articulation_variable_ids)
    isolated = set(residual.isolated_variable_ids)

    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="26" text-anchor="middle" '
        f'font-family="sans-serif" font-size="18">{escape(title)}</text>',
        f'<text x="{left_x}" y="46" text-anchor="middle" '
        f'font-family="sans-serif" font-size="12">residual legal decisions</text>',
        f'<text x="{right_x}" y="46" text-anchor="middle" '
        f'font-family="sans-serif" font-size="12">adjacent-layer crossing factors</text>',
    ]

    for factor_id in factor_ids:
        factor = factor_lookup[factor_id]
        for variable_id in factor.variable_ids:
            body.append(
                f'<line x1="{left_x}" y1="{variable_y[variable_id]}" '
                f'x2="{right_x}" y2="{factor_y[factor_id]}" '
                f'stroke="#888" stroke-opacity="0.28" stroke-width="1"/>'
            )

    for variable_id in variable_ids:
        y = variable_y[variable_id]
        is_articulation = variable_id in articulation
        is_isolated = variable_id in isolated
        stroke = "#b2182b" if is_articulation else "#222"
        fill = "#f3f4f6" if is_isolated else "white"
        stroke_width = 2.5 if is_articulation else 1.2
        label = variable_id.replace("order::", "order ").replace(
            "orient::", "orient "
        )
        body.append(
            f'<circle cx="{left_x}" cy="{y}" r="5" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{stroke_width}"/>'
        )
        body.append(
            f'<text x="{left_x - 10}" y="{y + 4}" text-anchor="end" '
            f'font-family="sans-serif" font-size="9">{escape(label)}</text>'
        )

    for factor_id in factor_ids:
        factor = factor_lookup[factor_id]
        y = factor_y[factor_id]
        label = f"{factor.species_left} ↔ {factor.species_right}"
        body.append(
            f'<rect x="{right_x - 4}" y="{y - 4}" width="8" height="8" '
            f'fill="white" stroke="#222" stroke-width="1.2"/>'
        )
        body.append(
            f'<text x="{right_x + 10}" y="{y + 4}" text-anchor="start" '
            f'font-family="sans-serif" font-size="9">{escape(label)}</text>'
        )

    body.append("</svg>")
    return "\n".join(body)


def _metric_table(rows: list[tuple[str, object]]) -> str:
    body = ["<table class='metrics'>"]
    for key, value in rows:
        body.append(
            "<tr><th>{}</th><td>{}</td></tr>".format(
                escape(str(key)),
                escape(str(value)),
            )
        )
    body.append("</table>")
    return "\n".join(body)


def render_validation_report_html(
    fixture: Fixture,
    result: ExactLayoutResult,
    *,
    solver: str,
    solver_details: dict[str, object] | None = None,
    hidden_native_state: LayoutState | None = None,
) -> str:
    """Create a self-contained visual/audit report for one solved fixture."""

    solver_details = {} if solver_details is None else solver_details
    audit = build_layout_audit(fixture, result)

    raw_graph = build_incidence_graph(fixture)
    raw_decomp = decompose_incidence_graph(raw_graph)
    structural = build_structural_projection(fixture)
    residual = build_residual_factorization(fixture)

    if result.optimality_status == "proven optimum":
        lower_bound = result.optimized_score.crossings
        upper_bound = result.optimized_score.crossings
    else:
        lower_bound = solver_details.get("lower_bound", "unknown")
        upper_bound = solver_details.get(
            "upper_bound", result.optimized_score.crossings
        )

    max_kernel_nodes = max(
        (
            len(kernel)
            for kernel in structural.decomposition.hard_kernels
        ),
        default=0,
    )

    metrics = _metric_table(
        [
            ("Fixture", fixture.fixture_id),
            ("Species", len(fixture.species_ids)),
            ("Chromosomes", len(fixture.chromosomes)),
            ("Homology groups", len(fixture.homology_ids)),
            ("Raw incidence cycle rank", sum(
                raw_graph.summarize_component(component).cycle_rank
                for component in raw_graph.connected_components()
            )),
            ("Structural bundles", len(structural.bundles)),
            ("Structural cycle rank", structural.cycle_rank_total),
            ("Structural hard kernels", len(
                structural.decomposition.hard_kernels
            )),
            ("Largest structural kernel (nodes)", max_kernel_nodes),
            ("Residual variables", len(residual.variables)),
            ("Residual crossing factors", len(residual.factors)),
            ("Residual objective components", residual.objective_component_count),
            (
                "Largest residual component (variables)",
                residual.max_objective_component_variables,
            ),
            (
                "Residual articulation variables",
                len(residual.articulation_variable_ids),
            ),
            (
                "Residual min-fill treewidth upper bound",
                residual.min_fill_treewidth_upper_bound,
            ),
            ("Initial crossings", result.initial_score.crossings),
            ("Normalized crossings", result.normalized_score.crossings),
            ("Optimized crossings", result.optimized_score.crossings),
            ("Lower bound", lower_bound),
            ("Upper bound", upper_bound),
            ("Optimality status", result.optimality_status),
            ("Solver", solver),
            ("States / nodes evaluated", result.states_evaluated),
            ("Input fingerprint", audit.input_fingerprint),
        ]
    )

    move_rows = []
    for move in audit.whole_chromosome_moves:
        move_rows.append(
            "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(
                escape(move.species_id),
                escape(move.chromosome_id),
                move.from_rank,
                move.to_rank,
            )
        )
    if not move_rows:
        move_rows.append(
            "<tr><td colspan='4'>No whole-chromosome rank changes.</td></tr>"
        )

    flip_text = (
        ", ".join(escape(item) for item in audit.whole_chromosome_flips)
        if audit.whole_chromosome_flips
        else "None"
    )

    native_section = ""
    if hidden_native_state is not None:
        native_crossings = score_crossings(
            fixture, hidden_native_state
        ).crossings
        native_svg = render_layout_state_svg(
            fixture,
            hidden_native_state,
            title="Hidden simulator-native baseline",
            crossing_count=native_crossings,
        )
        native_section = f"""
<section>
<h2>Simulator-only native baseline</h2>
<p class="note">
This panel is revealed only after solving. It is a validation baseline, not
solver input and not a claim that the simulator-native presentation is the
unique biologically correct layout.
</p>
<div class="svg-wrap">{native_svg}</div>
</section>
"""

    raw_svg = render_incidence_graph_svg(
        raw_graph,
        title="Raw chromosome ↔ homology incidence graph",
        core_nodes=raw_decomp.core_node_ids,
    )
    structural_svg = render_incidence_graph_svg(
        structural.graph,
        title="Chromosome-signature structural projection",
        core_nodes=structural.decomposition.core_node_ids,
    )
    comparison_svg = render_layout_comparison_svg(fixture, result)
    residual_svg = render_residual_factor_graph_svg(residual)

    details_json = escape(
        json.dumps(solver_details, indent=2, sort_keys=True)
    )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>SynTangle validation report — {escape(fixture.fixture_id)}</title>
<style>
body {{
  font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  margin: 2rem auto;
  max-width: 1500px;
  line-height: 1.4;
  color: #1d1d1f;
}}
h1, h2 {{ margin-top: 1.7rem; }}
.grid {{
  display: grid;
  grid-template-columns: minmax(320px, 0.8fr) minmax(520px, 1.2fr);
  gap: 1.5rem;
  align-items: start;
}}
.metrics, .moves {{
  border-collapse: collapse;
  width: 100%;
}}
.metrics th, .metrics td, .moves th, .moves td {{
  border-bottom: 1px solid #ddd;
  padding: 0.38rem 0.5rem;
  text-align: left;
  vertical-align: top;
}}
.metrics th {{ width: 48%; }}
.svg-wrap {{
  overflow-x: auto;
  border: 1px solid #ddd;
  background: white;
  padding: 0.4rem;
}}
.note {{
  padding: 0.7rem 0.9rem;
  border-left: 4px solid #777;
  background: #f5f5f5;
}}
pre {{
  white-space: pre-wrap;
  word-break: break-word;
  background: #f6f6f6;
  padding: 0.8rem;
  border: 1px solid #ddd;
}}
@media print {{
  body {{ max-width: none; margin: 0.5in; }}
  .svg-wrap {{ border: none; }}
}}
</style>
</head>
<body>
<h1>SynTangle validation report</h1>
<p>
This report audits the same extant biological input through decomposition,
legal layout optimization, and final visualization. Red graph nodes mark the
reported 2-core for that graph representation.
</p>

<div class="grid">
<section>
<h2>Run metrics</h2>
{metrics}
</section>
<section>
<h2>Transformations</h2>
<p><strong>Whole-chromosome flips:</strong> {flip_text}</p>
<table class="moves">
<tr><th>Species</th><th>Chromosome</th><th>Initial rank</th><th>Final rank</th></tr>
{''.join(move_rows)}
</table>
</section>
</div>

{native_section}

<section>
<h2>Tangled input → optimized layout</h2>
<div class="svg-wrap">{comparison_svg}</div>
</section>

<section>
<h2>Raw incidence structure</h2>
<p class="note">
This graph retains every homologous anchor. Its cycle rank can therefore
reflect repeated evidence density as well as chromosome-level structure.
</p>
<div class="svg-wrap">{raw_svg}</div>
</section>

<section>
<h2>Structural projection and hard core</h2>
<p class="note">
Homology groups with the same chromosome-incidence multiset are bundled only
for structural decomposition. Ordered block occurrences remain unchanged for
crossing and inversion evidence.
</p>
<div class="svg-wrap">{structural_svg}</div>
</section>

<section>
<h2>Residual decision / factor structure</h2>
<p class="note">
This is a derived computational graph, not a replacement for the biological
hypergraph. Left nodes are unresolved legal chromosome-order or GF(2)
orientation decisions; right nodes are adjacent-species crossing factors.
Red decision nodes are articulation variables in the residual primal graph.
Grey-filled decisions do not currently affect any crossing factor. Disconnected
residual pieces are exact candidates for the next solver-factorization layer.
</p>
<div class="svg-wrap">{residual_svg}</div>
</section>

<section>
<h2>Solver details</h2>
<pre>{details_json}</pre>
</section>
</body>
</html>
"""


def write_validation_report_html(
    fixture: Fixture,
    result: ExactLayoutResult,
    path: str,
    *,
    solver: str,
    solver_details: dict[str, object] | None = None,
    hidden_native_state: LayoutState | None = None,
) -> str:
    html = render_validation_report_html(
        fixture,
        result,
        solver=solver,
        solver_details=solver_details,
        hidden_native_state=hidden_native_state,
    )
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(html)
    return path
