import abc
import logging
from typing import Iterable

from gufe import AtomMapper, AtomMappingScorer
from gufe import LigandNetwork, Component

log = logging.getLogger(__name__)


# Todo: move to gufe Network planner

class NetworkPlanner(abc.ABC):

    def __init__(self, mapper: AtomMapper, scorer: AtomMappingScorer):
        """This class is an implementation for the LigandNetworkPlanner interface.
        It defines the std. class for a Konnektor LigandNetworkPlanner

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
        self.mapper = mapper
        self.scorer = scorer

    def __call__(self, *args, **kwargs) -> LigandNetwork:
        return self.generate_ligand_network(*args, **kwargs)

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
