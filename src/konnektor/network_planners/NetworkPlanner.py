# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import abc
from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, Component, LigandNetwork


def _validate_mappers(mappers) -> AtomMapper | Iterable[AtomMapper] | None:
    if isinstance(mappers, AtomMapper):
        return [mappers]
    elif isinstance(mappers, Iterable) and all(isinstance(m, AtomMapper) for m in mappers):
        return mappers
    elif mappers is None:
        return None
    else:
        raise ValueError("`mappers` must be of type `AtomMapper`.")


class NetworkPlanner(abc.ABC):
    def __init__(
        self,
        mappers: AtomMapper | Iterable[AtomMapper] | None,
        scorer: Callable[[AtomMapping], float] | None,
    ):
        """Abstract class for network planning, not to be called directly.

        Parameters
        ----------
        mappers : AtomMapper | Iterable[AtomMapper] | None
            AtomMapper(s) to use to propose mappings.
        scorer : Callable[[AtomMapping], float] | None
           Any callable which takes a AtomMapping and returns a float in [0,1].
        """

        self._mappers = _validate_mappers(mappers)
        self.scorer = scorer

    def __call__(self, *args, **kwargs) -> LigandNetwork:
        return self.generate_ligand_network(*args, **kwargs)

    @property
    def mappers(self) -> AtomMapper | Iterable[AtomMapper] | None:
        return self._mappers

    @mappers.setter
    def mappers(self, mappers: AtomMapper | Iterable[AtomMapper] | None):
        self._mappers = _validate_mappers(mappers)

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
        """Generate a LigandNetwork with the given Components as nodes and using this NetworkPlanner's mappers and scorer to create AtomMapping edges.

        Parameters
        ----------
        components : Iterable[Component]
            Component(s) to include in the Network

        Returns
        -------
        LigandNetwork
            LigandNetwork with Components as nodes and generated AtomMappings as edges.
        """
        raise NotImplementedError()
