# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
import logging
from collections.abc import Iterable

from gufe import AtomMapper, LigandNetwork

from ...network_planners._map_scoring import _score_mappings
from .._networkx_implementations import MstNetworkAlgorithm
from ._abstract_network_concatenator import NetworkConcatenator
from ...tests.network_planners.conftest import ligand_network_ab

log = logging.getLogger(__name__)


# TODO: check this algorithm again


class MstConcatenator(NetworkConcatenator):
    def __init__(
        self,
        mappers: AtomMapper | Iterable[AtomMapper] | None,
        scorer,
        n_connecting_edges: int = 2,
        n_processes: int = 1,
        _initial_edge_lister: NetworkConcatenator | None = None,  # TODO: remove this
    ):
        """
        A NetworkConcatenator that connects sub-networks with a minimum spanning tree.

        Parameters
        ----------
        mapper: AtomMapper
            AtomMapper(s) to use to propose mappings.
            If more than one AtomMapper is provided, the mapping with the lowest score (as scored by `scorer`) will be used.
        scorer: Callable[[AtomMapping], float] | None
            Callable which takes a AtomMapping and returns a float in [0,1].
        n_connecting_edges: int, optional
            Maximum number of edges to create when connecting the networks, by default 2.
        n_processes: int, optional
            Number of processes that can be used for the network generation, by default 1.
        """
        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=MstNetworkAlgorithm(),
            n_processes=n_processes,
            _initial_edge_lister=_initial_edge_lister,
        )
        self.n_connecting_edges = n_connecting_edges

    def _score_pair_edges(self, networkA: LigandNetwork, networkB: LigandNetwork) -> list:
        """Score every bipartite candidate edge between two sub-networks."""
        possible_edges = [(na, nb) for na in networkA.nodes for nb in networkB.nodes]
        return _score_mappings(
            possible_edges=possible_edges,
            scorer=self.scorer,
            mappers=self.mappers,
            n_processes=self.n_processes,
            show_progress=self.progress,
        )

    def _spanning_tree_pairs(
        self, pair_mappings: dict[tuple[int, int], list], n_networks: int
    ) -> list[tuple[int, int]]:
        """
        Build an MST over sub-networks to decide which sub-networks to join.

        Builds a graph whose nodes are the sub-networks, weighted by the best
        available score between each pair, and returns the pairs of a minimum
        spanning tree over the sub-networks.
        """
        # Create an "edge" between each sub-network pair
        super_edges = list(pair_mappings)
        # Weigh each sub-network connection by the score of the best possible
        # connection between those sub-networks.
        super_weights = [
            max(m.annotations["score"] for m in pair_mappings[e]) for e in
            super_edges
        ]
        # Create an MST where each node is a sub-network
        mst = self.network_generator.generate_network(
            super_edges, super_weights, n_edges=n_networks - 1
        )
        # Reorder the sub-network indices to match keys in pair_mappings
        subnetwork_pairs = [(min(e), max(e)) for e in mst.edges]
        return subnetwork_pairs

    def _select_connecting_edges(
            self, networkA: LigandNetwork, networkB: LigandNetwork, mappings: list
    ) -> list:
        """Pick up to `n_connecting_edges` mappings between two sub-networks."""
        ligands = list(networkA.nodes | networkB.nodes)
        edge_map = {
            (ligands.index(m.componentA), ligands.index(m.componentB)): m for m
            in mappings
        }
        edges = list(edge_map.keys())
        weights = [edge_map[k].annotations["score"] for k in edges]

        mg = self.network_generator.generate_network(
            edges, weights, n_edges=self.n_connecting_edges
        )
        return [
            edge_map[k] if k in edge_map else edge_map[tuple(list(k)[::-1])]
            for k in mg.edges
        ]

    def concatenate_networks(self, ligand_networks: Iterable[LigandNetwork]) -> LigandNetwork:
        """
        Concatenate the given networks.

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            LigandNetworks to concatenate.

        Returns
        -------
        LigandNetwork
            The concatenated LigandNetwork.
        """

        ligand_networks = list(ligand_networks)

        log.info(
            f"Number of edges in individual networks:\n"
            f"{sum(len(s.edges) for s in ligand_networks)}/"
            f"{[len(s.edges) for s in ligand_networks]}"
        )

        selected_edges = []
        selected_nodes = set()

        if len(ligand_networks) > 1:
            # Score candidate connecting edges for every pair of sub-networks
            pair_mappings = {}
            for i, j in itertools.combinations(range(len(ligand_networks)), 2):
                mappings = self._score_pair_edges(ligand_networks[i], ligand_networks[j])
                if mappings:
                    pair_mappings[(i, j)] = mappings

            # Identify which sub-networks to connect (MST over sub-networks)
            subnetwork_pairs = self._spanning_tree_pairs(pair_mappings, len(ligand_networks))

            # Connect each subnetwork pair with up to n_connecting_edges
            for i, j in subnetwork_pairs:
                connecting = self._select_connecting_edges(
                    ligand_networks[i],
                    ligand_networks[j],
                    pair_mappings[(i, j)]
                )
                log.info(f"Adding ConnectingEdges: {len(connecting)}")
                selected_edges.extend(connecting)

        # Constructed final Edges:
        # Add all old network edges:
        for network in ligand_networks:
            selected_edges.extend(network.edges)
            selected_nodes.extend(network.nodes)

        concat_network = LigandNetwork(edges=selected_edges, nodes=set(selected_nodes))
        log.info(f"Total Concatenated Edges: {len(selected_edges)}")

        if not concat_network.is_connected():
            raise RuntimeError("could not build a connected network!")

        return concat_network
