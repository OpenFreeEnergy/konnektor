# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import logging
from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, LigandNetwork

from .._networkx_implementations import MstNetworkAlgorithm
from ._abstract_network_concatenator import NetworkConcatenator
from .max_concatenator import MaxConcatenator

log = logging.getLogger(__name__)


# TODO: check this algorithm again
class CyclicConcatenator(NetworkConcatenator):
    def __init__(
        self,
        mappers: AtomMapper | Iterable[AtomMapper] | None,
        scorer: Callable[[AtomMapping], float] | None,
        n_connecting_cycles: int = 2,
        cycle_sizes: int | list[int] = 3,
        n_processes: int = 1,
        _initial_edge_lister: NetworkConcatenator | None = None,
    ):
        """

        Parameters
        ----------
        mappers: AtomMapper | Iterable[AtomMapper] | None
            AtomMapper(s) to use to propose mappings.
            If more than one AtomMapper is provided, the mapping with the lowest score (as scored by `scorer`) will be used.
        scorer: Callable[[AtomMapping], float] | None
            Callable which takes a AtomMapping and returns a float in [0,1].
        n_connecting_cycles: int, optional
            Minimum number of cycles to build between the networks, by default 2.
        cycle_sizes: Union[int, list[int]], optional
            Size of the cycles to build, can be an int or range of ints, by default 3.
        n_processes: int, optional
            Number of processes that can be used for the network generation, by default 1.
        _initial_edge_lister: NetworkConcatenator | None, optional
            The NetworkConcatenator to use if the NetworkConcatenator requires an initial set of edges, by default a MaxConcatenator with the provided `mappers`, `scorer` and `n_processes` will be used.
        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaxConcatenator(
                mappers=mappers, scorer=scorer, n_processes=n_processes
            )

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=MstNetworkAlgorithm(),
            n_processes=n_processes,
        )
        self.n_connecting_edges = n_connecting_cycles
        self.cycle_sizes = cycle_sizes

    def concatenate_networks(self, ligand_networks: Iterable[LigandNetwork]) -> LigandNetwork:
        """

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            an iterable of ligand networks, that shall be connected.
        n_connecting_edges: int
            number of edges, to connect the networks

        Returns
        -------
        LigandNetwork
            returns a concatenated LigandNetwork object,
            containing all networks.

        """
        raise NotImplementedError()
        # TODO: implement.

        # selected_edges = []
        # selected_nodes = []
        # for ligandNetworkA, ligandNetworkB in itertools.combinations(ligand_networks, 2):
        #     # Generate fully connected Bipartite Graph
        #     ligands = ligandNetworkA.nodes | ligandNetworkB.nodes
        #     fully_connected_graph = self._initial_edge_lister([ligandNetworkA, ligandNetworkB])
        #     bipartite_graph_mappings = list(fully_connected_graph.edges)

        #     # TODO Cycle Selection

        #     selected_edges.extend(selected_mappings)

        # # Constructed final Edges:
        # # Add all old network edges:
        # for network in ligand_networks:
        #     selected_edges.extend(network.edges)
        #     selected_nodes.extend(network.nodes)

        # concat_LigandNetwork = LigandNetwork(edges=selected_edges, nodes=set(selected_nodes))

        # log.info(f"Total Concatenated Edges: {len(selected_edges)}")

        # return concat_LigandNetwork
