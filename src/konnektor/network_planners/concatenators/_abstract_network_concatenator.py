# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import abc
import logging
from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, LigandNetwork

from .._networkx_implementations._abstract_network_algorithm import (
    _AbstractNetworkAlgorithm,
)
from ..NetworkPlanner import NetworkPlanner

log = logging.getLogger(__name__)


class NetworkConcatenator(NetworkPlanner):
    progress: bool = False
    n_processes: int

    def __init__(
        self,
        mappers: AtomMapper | Iterable[AtomMapper] | None,
        scorer: Callable[[AtomMapping], float] | None,
        network_generator: _AbstractNetworkAlgorithm | None,
        n_processes: int = 1,
        _initial_edge_lister=None,
    ):
        """Abstract class for network concatenation, not to be called directly.

        Parameters
        ----------
        mappers : AtomMapper | Iterable[AtomMapper] | None
            AtomMapper(s) to use to propose mappings.
            If more than one AtomMapper is provided, the mapping with the lowest score (as scored by `scorer`) will be used.
        scorer : Callable[[AtomMapping], float] | None
            Callable which takes a AtomMapping and returns a float in [0,1].
        n_processes: int, optional
            Number of processes that can be used for the network generation, by default 1.
        _initial_edge_lister: NetworkConcatenator | None, optional
            The NetworkConcatenator to use if the NetworkConcatenator requires an initial set of edges, by default None.
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

    def __call__(self, *args, **kwargs) -> LigandNetwork:
        return self.concatenate_networks(*args, **kwargs)

    @abc.abstractmethod
    def concatenate_networks(self, ligand_networks: Iterable[LigandNetwork]) -> LigandNetwork:
        """Concatenate the `ligand_networks` into a single LigandNetwork object.

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            LigandNetworks to concatenate.

        Returns
        -------
        LigandNetwork
            The concatenated LigandNetwork.
        """
        raise NotImplementedError()
