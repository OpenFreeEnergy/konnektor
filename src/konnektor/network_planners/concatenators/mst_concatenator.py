# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
import logging
from collections.abc import Iterable

from gufe import AtomMapper, LigandNetwork

from ...network_planners._map_scoring import _parallel_map_scoring
from .._networkx_implementations import MstNetworkAlgorithm
from ._abstract_network_concatenator import NetworkConcatenator

log = logging.getLogger(__name__)


# Todo: check this algorithm again


class MstConcatenator(NetworkConcatenator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
        scorer,
        n_connecting_edges: int = 2,
        n_processes: int = 1,
        _initial_edge_lister: NetworkConcatenator = None,
    ):
        """
        A NetworkConcatenator that connects two Networks with a Kruskal-like
        approach, up to the number of connecting edges.

        Parameters
        ----------
        mapper: AtomMapper
            the atom mapper is required, to define the connection between
            two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score
            between [0,1].
        n_connecting_edges: int, optional
            maximum number of connecting edges. (default: 2)
        n_processes: int
            number of processes that can be used for the network generation.
            (default: 1)
        """
        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=MstNetworkAlgorithm(),
            n_processes=n_processes,
            _initial_edge_lister=None,  # TODO: should this be _initial_edge_lister?
        )
        self.n_connecting_edges = n_connecting_edges

    def concatenate_networks(self, ligand_networks: Iterable[LigandNetwork]) -> LigandNetwork:
        """
        Concatenate the given networks.

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            an iterable of ligand networks, that shall be connected.

        Returns
        -------
        LigandNetwork
            returns a concatenated LigandNetwork object, containing all networks.

        """

        log.info(
            f"Number of edges in individual networks:\n"
            f"{sum([len(s.edges) for s in ligand_networks])}/"
            f"{[len(s.edges) for s in ligand_networks]}"
        )

        selected_edges = []
        selected_nodes = []
        for ligandNetworkA, ligandNetworkB in itertools.combinations(ligand_networks, 2):
            # Generate fully connected Bipartite Graph
            ligands = list(ligandNetworkA.nodes | ligandNetworkB.nodes)
            nodesA = ligandNetworkA.nodes
            nodesB = ligandNetworkB.nodes
            pedges = [(na, nb) for na in nodesA for nb in nodesB]

            bipartite_graph_mappings = _parallel_map_scoring(
                possible_edges=pedges,
                scorer=self.scorer,
                mappers=self.mappers,
                n_processes=self.n_processes,
                show_progress=self.progress,
            )

            # Find MST subset for Bipartite
            edge_map = {
                (ligands.index(m.componentA), ligands.index(m.componentB)): m
                for m in bipartite_graph_mappings
            }
            edges = list(edge_map.keys())
            weights = [edge_map[k].annotations["score"] for k in edges]

            mg = self.network_generator.generate_network(
                edges, weights, n_edges=self.n_connecting_edges
            )

            selected_mappings = [
                edge_map[k] if (k in edge_map) else edge_map[tuple(list(k)[::-1])] for k in mg.edges
            ]

            log.info(f"Adding ConnectingEdges: {len(selected_mappings)}")

            # Add network connecting edges
            selected_edges.extend(selected_mappings)

        # Constructed final Edges:
        # Add all old network edges:
        for network in ligand_networks:
            selected_edges.extend(network.edges)
            selected_nodes.extend(network.nodes)

        concat_LigandNetwork = LigandNetwork(edges=selected_edges, nodes=set(selected_nodes))
        log.info(f"Total Concatenated Edges: {len(selected_edges)}")

        return concat_LigandNetwork
