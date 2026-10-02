from __future__ import annotations

from html import escape

from .layout import ExactLayoutResult, LayoutState
from .model import Fixture


def _panel_geometry(
    fixture: Fixture,
    state: LayoutState,
    panel_x: float,
    panel_width: float,
    top: float,
    row_gap: float,
) -> tuple[dict[str, float], dict[object, tuple[float, float, float]]]:
    species_y = {
        species: top + index * row_gap
        for index, species in enumerate(fixture.species_ids)
    }
    chrom_geom: dict[object, tuple[float, float, float]] = {}

    usable_width = panel_width - 100.0
    for species in fixture.species_ids:
        refs = state.chromosome_order[species]
        lengths = {
            chrom.ref: chrom.length
            for chrom in fixture.chromosomes
            if chrom.ref.species_id == species
        }
        if not refs:
            continue
        gap = 8.0
        total_gap = gap * max(0, len(refs) - 1)
        total_length = sum(lengths[ref] for ref in refs)
        scale = (usable_width - total_gap) / total_length if total_length else 1.0
        x = panel_x + 75.0
        for ref in refs:
            width = lengths[ref] * scale
            chrom_geom[ref] = (x, width, species_y[species])
            x += width + gap

    return species_y, chrom_geom


def _anchor_x(fixture: Fixture, state: LayoutState, chrom_geom, chromosome, block) -> float:
    x, width, _ = chrom_geom[chromosome.ref]
    midpoint = (block.start + block.end) / 2.0
    fraction = midpoint / chromosome.length
    if state.chromosome_orientation[chromosome.ref] == -1:
        fraction = 1.0 - fraction
    return x + fraction * width


def _unambiguous_links(fixture: Fixture, species1: str, species2: str):
    by_species: dict[str, dict[str, list[tuple[object, object]]]] = {}
    for chromosome in fixture.chromosomes:
        for block in chromosome.blocks:
            by_species.setdefault(chromosome.ref.species_id, {}).setdefault(
                block.homology_id, []
            ).append((chromosome, block))

    homologies = sorted(
        set(by_species.get(species1, {}))
        & set(by_species.get(species2, {}))
    )
    links = []
    for homology_id in homologies:
        left = by_species[species1][homology_id]
        right = by_species[species2][homology_id]
        if len(left) == 1 and len(right) == 1:
            links.append((homology_id, left[0], right[0]))
    return links


def _render_panel(
    fixture: Fixture,
    state: LayoutState,
    panel_x: float,
    panel_width: float,
    top: float,
    row_gap: float,
    title: str,
    crossing_count: int,
) -> list[str]:
    species_y, chrom_geom = _panel_geometry(
        fixture, state, panel_x, panel_width, top, row_gap
    )
    svg: list[str] = []
    svg.append(
        f'<text x="{panel_x + panel_width / 2:.1f}" y="28" '
        f'text-anchor="middle" font-size="18" font-family="sans-serif">'
        f'{escape(title)} — crossings: {crossing_count}</text>'
    )

    for species1, species2 in zip(fixture.species_ids, fixture.species_ids[1:]):
        for _, (chrom1, block1), (chrom2, block2) in _unambiguous_links(
            fixture, species1, species2
        ):
            x1 = _anchor_x(fixture, state, chrom_geom, chrom1, block1)
            x2 = _anchor_x(fixture, state, chrom_geom, chrom2, block2)
            y1 = species_y[species1]
            y2 = species_y[species2]
            svg.append(
                f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
                f'stroke="#777" stroke-opacity="0.25" stroke-width="1"/>'
            )

    for species in fixture.species_ids:
        y = species_y[species]
        svg.append(
            f'<text x="{panel_x + 8:.1f}" y="{y + 5:.1f}" '
            f'font-size="13" font-family="sans-serif">{escape(species)}</text>'
        )
        for ref in state.chromosome_order[species]:
            x, width, _ = chrom_geom[ref]
            orientation = state.chromosome_orientation[ref]
            svg.append(
                f'<rect x="{x:.2f}" y="{y - 8:.2f}" width="{max(width, 1):.2f}" '
                f'height="16" rx="3" fill="white" stroke="black" stroke-width="1.4"/>'
            )
            label = ref.chromosome_id if orientation == 1 else f"{ref.chromosome_id} ↔"
            svg.append(
                f'<text x="{x + width / 2:.2f}" y="{y - 13:.2f}" text-anchor="middle" '
                f'font-size="9" font-family="sans-serif">{escape(label)}</text>'
            )
    return svg


def render_layout_comparison_svg(
    fixture: Fixture,
    result: ExactLayoutResult,
    *,
    panel_width: int = 760,
    row_gap: int = 120,
) -> str:
    """Render initial and optimized layouts from the exact same biological input."""

    top = 85
    height = top + row_gap * max(1, len(fixture.species_ids) - 1) + 80
    width = panel_width * 2

    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]
    body.extend(
        _render_panel(
            fixture,
            result.initial_state,
            0,
            panel_width,
            top,
            row_gap,
            "Initial",
            result.initial_score.crossings,
        )
    )
    body.append(
        f'<line x1="{panel_width}" y1="40" x2="{panel_width}" y2="{height - 25}" '
        f'stroke="#ccc" stroke-width="1"/>'
    )
    body.extend(
        _render_panel(
            fixture,
            result.optimized_state,
            panel_width,
            panel_width,
            top,
            row_gap,
            "Optimized",
            result.optimized_score.crossings,
        )
    )
    body.append("</svg>")
    return "\n".join(body)


def write_layout_comparison_svg(
    fixture: Fixture,
    result: ExactLayoutResult,
    path: str,
) -> str:
    svg = render_layout_comparison_svg(fixture, result)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(svg)
    return path
