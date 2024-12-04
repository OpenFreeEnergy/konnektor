# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import abc
import logging
from typing import Iterable, Union

from gufe import AtomMapper
from gufe import LigandNetwork

from ..NetworkPlanner import NetworkPlanner
from .._networkx_implementations._abstract_network_algorithm import (
    _AbstractNetworkAlgorithm,
)

log = logging.getLogger(__name__)


class NetworkConcatenator(NetworkPlanner):
    progress: bool = False
    n_processes: int

    def __init__(
        self,
        mappers: Union[AtomMapper, Iterable[AtomMapper]],
        scorer,
        network_generator: _AbstractNetworkAlgorithm,
        n_processes: int = 1,
        _initial_edge_lister=None,
    ):
        """Base Class for the NetworkConcatenator classes.
         It defines the std. class for a Konnektor NetworkConcatenator.

        Parameters
        ----------
        mappers : AtomMapper
            the AtomMappers to use to propose mappings.  At least 1 required,
            but many can be given, in which case all will be tried to find the
            lowest score edges
        scorer : AtomMappingScorer
            any callable which takes a AtomMapping and returns a float
        n_processes: int, optional
            number of processes that can be used for the network generation.
            (default: 1)
        _initial_edge_lister: LigandNetworkPlanner, optional
            this LigandNetworkPlanner is used to give the initial set of edges.
             For standard usage, the Maximal NetworPlanner is used.
            However in large scale approaches, it might be interesting to use
             the heuristicMaximalNetworkPlanner. (default: None)

        """

        # generic Network_Planner attribsd
        super().__init__(mappers=mappers, scorer=scorer)

        # Konnektor specific variables
        self.network_generator = network_generator
        self.n_processes = n_processes
        self._initial_edge_lister = _initial_edge_lister

        # pass on the parallelization to the edge lister
        # edge lister performs usually the most expensive task!
        # So parallelization is most important here.
        if self._initial_edge_lister is not None and hasattr(
            self._initial_edge_lister, "nprocesses"
        ):
            self.n_processes = n_processes

    def __call__(self, *args, **kwargs) -> LigandNetwork:
        return self.concatenate_networks(*args, **kwargs)

    @abc.abstractmethod
    def concatenate_networks(
        self, ligand_networks: Iterable[LigandNetwork]
    ) -> LigandNetwork:
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
