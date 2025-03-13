# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
import logging
from collections.abc import Iterable

from gufe import AtomMapper, LigandNetwork

from ...network_planners._map_scoring import _parallel_map_scoring, _serial_map_scoring
from ._abstract_network_concatenator import NetworkConcatenator

log = logging.getLogger(__name__)


class MaxConcatenator(NetworkConcatenator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
        scorer,
        n_processes: int = 1,
        show_progress: bool = False,
    ):
        """
        A NetworkConcatenator that connects two Networks with all possible
        mappings. This is usually most useful for initial edge listing.

        Parameters
        ----------
        mappers: AtomMapper
            the atom mapper is required to define the connection
            between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a
            score between [0,1].
        n_processes: int
            number of processes that can be used for the network generation.
            (default: 1)
        show_progress: bool
            show progress bar
        """

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=None,
            n_processes=n_processes,
        )
        self.progress = show_progress

    def concatenate_networks(self, ligand_networks: Iterable[LigandNetwork]) -> LigandNetwork:
        """

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            An iterable of LigandNetworks to connect.

        Returns
        -------
        LigandNetwork
            returns a concatenated LigandNetwork object, containing all
             networks and all possible edges, connecting them.

        """

        log.info(
            f"Number of edges in individual networks:\n"
            f"{sum([len(s.edges) for s in ligand_networks])}/"
            f"{[len(s.edges) for s in ligand_networks]}"
        )

        selected_edges = []
        selected_nodes = []
        for ligandNetworkA, ligandNetworkB in itertools.combinations(ligand_networks, 2):
            # Generate Full Bipartite Graph
            nodesA = ligandNetworkA.nodes
            nodesB = ligandNetworkB.nodes
            pedges = [(na, nb) for na in nodesA for nb in nodesB]

            if self.n_processes > 1:
                bipartite_graph_mappings = _parallel_map_scoring(
                    possible_edges=pedges,
                    scorer=self.scorer,
                    mappers=self.mappers,
                    n_processes=self.n_processes,
                    show_progress=self.progress,
                )

            else:  # serial variant
                bipartite_graph_mappings = _serial_map_scoring(
                    possible_edges=pedges,
                    scorer=self.scorer,
                    mappers=self.mappers,
                    n_edges_to_score=len(pedges),
                    show_progress=self.progress,
                )

            # Add network connecting edges
            selected_edges.extend(bipartite_graph_mappings)

        # Constructed final Edges:
        # Add all old network edges:
        for network in ligand_networks:
            selected_edges.extend(network.edges)
            selected_nodes.extend(network.nodes)

        concat_LigandNetwork = LigandNetwork(edges=selected_edges, nodes=set(selected_nodes))

        log.info(f"Total Concatenated Edges: {len(selected_edges)} ")

        if not concat_LigandNetwork.is_connected():
            raise RuntimeError("could not build a connected network!")

        return concat_LigandNetwork
