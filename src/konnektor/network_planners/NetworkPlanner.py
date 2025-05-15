# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import abc
import logging
from collections.abc import Iterable

from gufe import AtomMapper, Component, LigandNetwork

log = logging.getLogger(__name__)


# Todo: move to gufe Network planner


class NetworkPlanner(abc.ABC):
    def __init__(self, mappers: AtomMapper | list[AtomMapper], scorer):
        """This class is an implementation for the NetworkPlanner interface.
        It defines the std. class for a Konnektor NetworkPlanner

        Parameters
        ----------
        mapper : AtomMapper
            the AtomMappers to use to propose mappings.  At least 1 required,
            but many can be given, in which case all will be tried to find the
            lowest score edges
        scorer : AtomMappingScorer
            any callable which takes a AtomMapping and returns a float
        """

        # generic Network_Planner attribs
        if isinstance(mappers, AtomMapper):
            self._mappers = [mappers]
        elif isinstance(mappers, Iterable) and all(isinstance(m, AtomMapper) for m in mappers):
            self._mappers = mappers
        elif mappers is None:
            self._mappers = None
        else:
            raise ValueError("Atom mappers are not the required type!")

        self.scorer = scorer

    def __call__(self, *args, **kwargs) -> LigandNetwork:
        return self.generate_ligand_network(*args, **kwargs)

    @property
    def mappers(self) -> list[AtomMapper]:
        return self._mappers

    @mappers.setter
    def mappers(self, mappers: AtomMapper | list[AtomMapper]):
        if mappers is AtomMapper:
            self._mappers = [mappers]
        elif isinstance(mappers, Iterable) and all(isinstance(m, AtomMapper) for m in mappers):
            self._mappers = mappers

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
        """Plan a Network which connects all ligands following a given algorithm cost

        Parameters
        ----------
        components : Iterable[Component]
        the ligands to include in the Network

        Returns
        -------
        LigandNetwork
            the resulting ligand network.
        """
        raise NotImplementedError()
