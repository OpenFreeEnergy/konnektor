import itertools
import logging
from typing import Iterable

from gufe import AtomMapper, LigandNetwork

from ._abstract_network_concatenator import NetworkConcatenator
from ..generators._parallel_mapping_pattern import _parallel_map_scoring

log = logging.getLogger(__name__)


class MaxConcatenate(NetworkConcatenator):
    def __init__(self, mapper: AtomMapper, scorer, nprocesses: int = 1,
                 show_progress: bool = False):
        """
        This concatenator is connnecting two Networks with all possible
         mappings. This is usually most useful for initial edge listing.

        Parameters
        ----------
        mapper: AtomMapper
            the atom mapper is required, to define the connection
             between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a
            score between [0,1].
        n_connecting_edges: int, optional
            number of connecting edges. (default: 3)
        nprocesses: int
            number of processes that can be used for the network generation.
            (default: 1)
        show_progress: bool
            show progress bar
        """

        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=None,
                         nprocesses=nprocesses)
        self.progress = show_progress

    def concatenate_networks(self, ligand_networks: Iterable[
        LigandNetwork]) -> LigandNetwork:
        """

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            an iterable of ligand networks, that shall be connected.

        Returns
        -------
        LigandNetwork
            returns a concatenated LigandNetwork object, containing all
             networks and all possible edges, connecting them.

        """

        log.info("Number of edges in individual networks:\n " + str(
            sum([len(s.edges) for s in ligand_networks])) +
                 "/ " + str([len(s.edges) for s in ligand_networks]))
        selected_edges = []
        selected_nodes = []
        for ligandNetworkA, ligandNetworkB in itertools.combinations(
                ligand_networks, 2):
            # Generate Full Bipartite Graph
            nodesA = ligandNetworkA.nodes
            nodesB = ligandNetworkB.nodes
            pedges = [(na, nb) for na in nodesA for nb in nodesB]

            bipartite_graph_mappings = _parallel_map_scoring(
                possible_edges=pedges,
                scorer=self.scorer,
                mapper=self.mapper, n_processes=self.nprocesses,
                show_progress=self.progress)

            # Add network connecting edges
            selected_edges.extend(bipartite_graph_mappings)

        # Constructed final Edges:
        # Add all old network edges:
        for network in ligand_networks:
            selected_edges.extend(network.edges)
            selected_nodes.extend(network.nodes)

        concat_LigandNetwork = LigandNetwork(edges=selected_edges,
                                             nodes=set(selected_nodes))

        log.info("Total Concatenated Edges:  " + str(len(selected_edges)))

        return concat_LigandNetwork
