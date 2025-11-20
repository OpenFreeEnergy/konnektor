from collections.abc import Callable
from typing import TypeAlias

from gufe import AtomMapping

# TODO: this will be implemented in gufe
AtomMappingScorer: TypeAlias = Callable[[AtomMapping], float]
"""A function that takes in an ``AtomMapping`` and returns a float between 0 (worst) and 1 (best) indicating the quality of the mapping."""
