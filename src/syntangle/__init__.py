"""SynTangle core package."""

from .cycles import CycleBasisElement, fundamental_cycle_basis
from .decomposition import GraphDecomposition, decompose_incidence_graph
from .fixtures import fixture_from_dict, load_fixture
from .incidence import analyze_fixture, build_incidence_graph
from .model import Fixture, FixtureValidationError
from .ordering import OrderingConstraintState, derive_ordering_constraints
from .orientation import OrientationResult, solve_orientation_constraints

__all__ = [
    "CycleBasisElement",
    "Fixture",
    "FixtureValidationError",
    "GraphDecomposition",
    "OrderingConstraintState",
    "OrientationResult",
    "analyze_fixture",
    "build_incidence_graph",
    "decompose_incidence_graph",
    "derive_ordering_constraints",
    "fixture_from_dict",
    "fundamental_cycle_basis",
    "load_fixture",
    "solve_orientation_constraints",
]
