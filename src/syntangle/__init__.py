"""SynTangle core package."""

from .fixtures import fixture_from_dict, load_fixture
from .incidence import analyze_fixture, build_incidence_graph
from .model import Fixture, FixtureValidationError

__all__ = [
    "Fixture",
    "FixtureValidationError",
    "analyze_fixture",
    "build_incidence_graph",
    "fixture_from_dict",
    "load_fixture",
]
