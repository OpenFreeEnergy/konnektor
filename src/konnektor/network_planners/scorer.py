from collections.abc import Callable
from typing import TypeAlias

from gufe import AtomMapping

AtomMappingScorer: TypeAlias = Callable[
    [AtomMapping], float
]  # TODO: this will be implemented in gufe
"""A function that takes in an ``AtomMapping`` and returns a float between 0 and 1."""
