"""SynTangle core package."""

from .decomposition import GraphDecomposition, decompose_incidence_graph
from .fixtures import fixture_from_dict, load_fixture
from .incidence import analyze_fixture, build_incidence_graph
from .model import Fixture, FixtureValidationError
from .orientation import OrientationResult, solve_orientation_constraints

__all__ = [
    "Fixture",
    "FixtureValidationError",
    "GraphDecomposition",
    "OrientationResult",
    "analyze_fixture",
    "build_incidence_graph",
    "decompose_incidence_graph",
    "fixture_from_dict",
    "load_fixture",
    "solve_orientation_constraints",
]
