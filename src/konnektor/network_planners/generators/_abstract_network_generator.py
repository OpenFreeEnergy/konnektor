# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import abc
import logging
from typing import Callable, Iterable, Union

from gufe import AtomMapper, AtomMapping, Component, LigandNetwork

from konnektor.network_planners._networkx_implementations import (
    _AbstractNetworkAlgorithm,
)

from ..NetworkPlanner import NetworkPlanner

log = logging.getLogger(__name__)


class NetworkGenerator(NetworkPlanner):
    progress: bool = False
    n_processes: int

    def __init__(
        self,
        mappers: Union[AtomMapper, list[AtomMapper]],
        scorer: Callable[[AtomMapping], float],
        network_generator: _AbstractNetworkAlgorithm,  # TODO: rename this to network_algorithm?
        n_processes: int = 1,
        progress: bool = False,
        _initial_edge_lister=None,
    ):
        """This class is an implementation for the NetworkGenerator interface.
        It defines the std. class for a Konnektor NetworkGenerator.

        Parameters
        ----------
        mappers : AtomMapper
            the AtomMappers to use to propose mappings.  At least 1 required,
            but many can be given, in which case all will be tried to find the
            lowest score edges
        scorer : AtomMappingScorer
            Any callable which takes a AtomMapping and returns a float
        network_generator: the network algorithm to use
        n_processes: int, optional
            Number of processes that can be used for the network generation. (default: 1)
        progress: bool, optional
            If `True`, displays a progress bar. (default: False)
        _initial_edge_lister: NetworkPlanner, optional
            The NetworkPlanner to use to create the initial set of edges. For standard usage, the MaximalNetworkPlanner is used.
            However in large scale approaches, it might be interesting to use the heuristicMaximalNetworkPlanner. (default: None)

        """
        # generic NetworkPlanner attribs
        super().__init__(mappers=mappers, scorer=scorer)

        # Konnektor specific variables
        self.network_generator = network_generator
        self.n_processes = n_processes

        self._initial_edge_lister = _initial_edge_lister

        # pass on the parallelization to the edge lister
        # edge lister performs usually the most expensive task!
        # So parallelization is most important here.
        if self._initial_edge_lister is not None and hasattr(
            self._initial_edge_lister, "n_processes"
        ):
            self.n_processes = n_processes
        if self._initial_edge_lister is not None and hasattr(self._initial_edge_lister, "progress"):
            self._initial_edge_lister._progress = progress
        self._progress = progress

    @property
    def progress(self) -> bool:
        """
        shows a progress bar if True
        """
        return self._progress

    @progress.setter
    def progress(self, progress: bool):
        self._progress = progress
        if self._initial_edge_lister is not None and hasattr(self._initial_edge_lister, "progress"):
            self._initial_edge_lister._progress = progress

    @abc.abstractmethod
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
