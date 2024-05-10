import abc
import logging
from typing import Iterable

from gufe import AtomMapper
from gufe import LigandNetwork

from ..NetworkPlanner import NetworkPlanner
from .._networkx_implementations._abstract_network_generator import \
    _AbstractNetworkGenerator

log = logging.getLogger(__name__)


class NetworkConcatenator(NetworkPlanner):
    progress: bool = False
    nprocesses: int

    def __init__(self, mapper: AtomMapper, scorer,
                 network_generator: _AbstractNetworkGenerator,
                 nprocesses: int = 1,
                 _initial_edge_lister=None):
        """This class is an implementation for the LigandNetworkPlanner
        interface. It defines the std. class for a
        Konnektor LigandNetworkPlanner

        Parameters
        ----------
        mapper : AtomMapper
            the AtomMappers to use to propose mappings.  At least 1 required,
            but many can be given, in which case all will be tried to find the
            lowest score edges
        scorer : AtomMappingScorer
            any callable which takes a AtomMapping and returns a float
        nprocesses: int, optional
            number of processes that can be used for the network generation.
            (default: 1)
        _initial_edge_lister: LigandNetworkPlanner, optional
            this LigandNetworkPlanner is used to give the initial set of edges.
             For standard usage, the Maximal NetworPlanner is used.
            However in large scale approaches, it might be interesting to use
             the heuristicMaximalNetworkPlanner. (default: None)

        """

        # generic Network_Planner attribsd
        super().__init__(mapper=mapper, scorer=scorer)

        # Konnektor specific variables
        self.network_generator = network_generator
        self.nprocesses = nprocesses
        self._initial_edge_lister = _initial_edge_lister

        # pass on the parallelization to the edge lister
        # edge lister performs usually the most expensive task!
        # So parallelization is most important here.
        if self._initial_edge_lister is not None and hasattr(
                self._initial_edge_lister, "nprocesses"):
            self.nprocesses = nprocesses

    def __call__(self, *args, **kwargs) -> LigandNetwork:
        return self.concatenate_networks(*args, **kwargs)

    @abc.abstractmethod
    def concatenate_networks(self, ligand_networks: Iterable[
        LigandNetwork]) -> LigandNetwork:
        """

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            an iterable of ligand networks, that shall be connected.

        Returns
        -------
        LigandNetwork
            returns a concatenated LigandNetwork object, containing all networks.

        """
