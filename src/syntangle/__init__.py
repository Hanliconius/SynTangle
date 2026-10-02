"""SynTangle core package."""

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

__all__ = [
    "CrossingScore",
    "CycleBasisElement",
    "ExactLayoutResult",
    "Fixture",
    "FixtureValidationError",
    "GraphDecomposition",
    "LayoutState",
    "OrderingConstraintState",
    "OrientationResult",
    "analyze_fixture",
    "build_incidence_graph",
    "canonicalize_component_order",
    "decompose_incidence_graph",
    "derive_ordering_constraints",
    "exact_optimize_small",
    "fixture_from_dict",
    "fundamental_cycle_basis",
    "initial_layout_state",
    "load_fixture",
    "score_crossings",
    "solve_orientation_constraints",
]
