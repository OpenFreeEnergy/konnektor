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
            self._initial_edge_lister, "nprocesses"
        ):
            self.n_processes = n_processes

    def __call__(self, *args, **kwargs) -> LigandNetwork:
        return self.concatenate_networks(*args, **kwargs)

    @staticmethod
    def _normalize_excluded_edges(
        exclude_edges: Iterable[AtomMapping] | None,
    ) -> set[frozenset]:
        """Normalize `exclude_edges` into a set of undirected component-pairs."""
        return {
            frozenset((mapping.componentA, mapping.componentB)) for mapping in (exclude_edges or ())
        }

    @staticmethod
    def _generate_bipartite_edges(
        networkA: LigandNetwork,
        networkB: LigandNetwork,
        exclude: set[frozenset],
    ) -> list[tuple]:
        """Generate allowed edges between two ligand networks."""
        return [
            (ligandA, ligandB)
            for ligandA in networkA.nodes
            for ligandB in networkB.nodes
            if frozenset((ligandA, ligandB)) not in exclude
        ]

    @abc.abstractmethod
    def concatenate_networks(
        self,
        ligand_networks: Iterable[LigandNetwork],
        exclude_edges: Iterable[AtomMapping] | None = None,
    ) -> LigandNetwork:
        """Concatenate the `ligand_networks` into a single LigandNetwork object.

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            LigandNetworks to concatenate.
        exclude_edges: Iterable[AtomMapping], optional
            Mappings that cannot be proposed as new connections which is useful for excluding edges
            that had already failed. If excluding these edges leaves the network unbridgeable, an error is raised.
            Default: None

        Returns
        -------
        LigandNetwork
            The concatenated LigandNetwork.
        """
        raise NotImplementedError()
