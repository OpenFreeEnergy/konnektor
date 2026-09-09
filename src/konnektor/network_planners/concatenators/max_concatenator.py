# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
import logging
from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, LigandNetwork

from ...network_planners._map_scoring import _score_mappings
from ._abstract_network_concatenator import NetworkConcatenator

log = logging.getLogger(__name__)


class MaxConcatenator(NetworkConcatenator):
    def __init__(
        self,
        mappers: AtomMapper | Iterable[AtomMapper] | None,
        scorer: Callable[[AtomMapping], float] | None,
        n_processes: int = 1,
        show_progress: bool = False,
    ):
        """
        A NetworkConcatenator that connects a set of LigandNetworks with all possible edges.
        This is usually most useful for initial edge listing.

        Parameters
        ----------
        mappers : AtomMapper | Iterable[AtomMapper] | None
            AtomMapper(s) to use to propose mappings. If more than one AtomMapper is provided, all will be tried to find the
            lowest score for each edges.
        scorer : Callable[[AtomMapping], float] | None
           Any callable which takes a AtomMapping and returns a float in [0,1].
        n_processes: int, optional
            Number of processes that can be used for the network generation, by default 1.
        show_progress:  bool, optional
            If True, a progress bar will be displayed, by default False.
        """

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=None,
            n_processes=n_processes,
        )
        self.progress = show_progress

    def _concatenate_networks(
        self,
        ligand_networks: Iterable[LigandNetwork],
        exclude: set[frozenset],
    ) -> LigandNetwork:
        """
        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            An iterable of LigandNetworks to connect.
        exclude : set[frozenset]
            Unordered ligand pairs that must not be proposed as new connections.

        Returns
        -------
        LigandNetwork
            The concatenated LigandNetwork with all possible nodes connected by edges.
        """
        log.info(
            f"Number of edges in individual networks:\n"
            f"{sum([len(s.edges) for s in ligand_networks])}/"
            f"{[len(s.edges) for s in ligand_networks]}"
        )

        selected_edges = []
        selected_nodes = set()
        for networkA, networkB in itertools.combinations(ligand_networks, 2):
            # Generate Full Bipartite Graph, excluding specified edges
            possible_edges = self._generate_bipartite_edges(networkA, networkB, exclude)
            bipartite_graph_mappings = _score_mappings(
                possible_edges=possible_edges,
                scorer=self.scorer,
                mappers=self.mappers,
                n_processes=self.n_processes,
                show_progress=self.progress,
            )
            # Add network connecting edges
            selected_edges.extend(bipartite_graph_mappings)

        # Constructed final Edges:
        # Add all old network edges:
        for network in ligand_networks:
            selected_edges.extend(network.edges)
            selected_nodes.update(network.nodes)

        concat_network = LigandNetwork(edges=selected_edges, nodes=selected_nodes)

        log.info(f"Total Concatenated Edges: {len(selected_edges)} ")

        return concat_network
