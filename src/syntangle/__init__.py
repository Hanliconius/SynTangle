"""SynTangle core package."""

from .audit import LayoutAudit, build_layout_audit, fixture_fingerprint
from .cycles import CycleBasisElement, fundamental_cycle_basis
from .decomposition import GraphDecomposition, decompose_incidence_graph
from .fixtures import fixture_from_dict, load_fixture
from .incidence import analyze_fixture, build_incidence_graph
from .layout import (
    CrossingScore,
    ExactLayoutResult,
    LayoutState,
    canonicalize_component_order,
    exact_optimize_small,
    initial_layout_state,
    score_crossings,
)
from .model import Fixture, FixtureValidationError
from .ordering import OrderingConstraintState, derive_ordering_constraints
from .orientation import OrientationResult, solve_orientation_constraints
from .visualize import render_layout_comparison_svg, write_layout_comparison_svg

__all__ = [
    "CrossingScore",
    "CycleBasisElement",
    "ExactLayoutResult",
    "Fixture",
    "FixtureValidationError",
    "GraphDecomposition",
    "LayoutAudit",
    "LayoutState",
    "OrderingConstraintState",
    "OrientationResult",
    "analyze_fixture",
    "build_incidence_graph",
    "build_layout_audit",
    "canonicalize_component_order",
    "decompose_incidence_graph",
    "derive_ordering_constraints",
    "exact_optimize_small",
    "fixture_fingerprint",
    "fixture_from_dict",
    "fundamental_cycle_basis",
    "initial_layout_state",
    "load_fixture",
    "render_layout_comparison_svg",
    "score_crossings",
    "solve_orientation_constraints",
    "write_layout_comparison_svg",
]
