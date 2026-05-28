# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import abc
from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, Component, LigandNetwork

from konnektor.network_planners._networkx_implementations import (
    _AbstractNetworkAlgorithm,
)

from ..NetworkPlanner import NetworkPlanner


class NetworkGenerator(NetworkPlanner):
    def __init__(
        self,
        mappers: AtomMapper | Iterable[AtomMapper] | None,
        scorer: Callable[[AtomMapping], float] | None,
        # TODO: rename this to network_algorithm and make base class not private?
        network_generator: _AbstractNetworkAlgorithm | None,
        n_processes: int = 1,
        progress: bool = False,
        _initial_edge_lister=None,  # TODO: why is this private, can we deprecate and rename?
    ):
        """Abstract class for network generation, not to be called directly.

        Parameters
        ----------
        mappers : AtomMapper | Iterable[AtomMapper] | None
            AtomMapper(s) to use to propose mappings. If more than one AtomMapper is provided, the mapping with the lowest score (as scored by `scorer`) will be used.
        scorer : Callable[[AtomMapping], float] | None
            Callable which takes a AtomMapping and returns a float in [0,1].
        network_generator : _AbstractNetworkAlgorithm | None.
            Algorithm to use when generating the network.
        n_processes : int, optional
            Number of processes that can be used for the network generation, by default 1.
        progress : bool, optional
            If True, display a progress bar, by default False.
        _initial_edge_lister : NetworkGenerator | None
            The NetworkGenerator to use if the NetworkGenerator requires an initial set of edges, by default None.
        """

        super().__init__(mappers=mappers, scorer=scorer)

        self.network_generator = network_generator
        self.n_processes = n_processes
        self._initial_edge_lister = _initial_edge_lister

        # pass on the parallelization to the edge lister
        # edge listing is usually the most expensive task,
        # so parallelization is important here.
        if self._initial_edge_lister is not None and hasattr(
            self._initial_edge_lister, "n_processes"
        ):
            self.n_processes = n_processes
        if self._initial_edge_lister is not None and hasattr(self._initial_edge_lister, "progress"):
            self._initial_edge_lister._progress = progress
        self._progress = progress

    @property
    def progress(self) -> bool:
        """If True, displays a progress bar when generating the LigandNetwork."""
        return self._progress

    @progress.setter
    def progress(self, progress: bool):
        self._progress = progress
        if self._initial_edge_lister is not None and hasattr(self._initial_edge_lister, "progress"):
            self._initial_edge_lister._progress = progress

    @abc.abstractmethod
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
