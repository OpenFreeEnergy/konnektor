# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import logging

import networkx as nx

from ._abstract_network_algorithm import _AbstractNetworkConcatenator

log = logging.getLogger(__name__)


class MatchingConcatenatAlgorithm(_AbstractNetworkConcatenator):

    def concatenate_networks(
        self,
        nodesA: list[int],
        nodesB: list[int],
        edges: list[tuple[int, int]],
        weights: list[float],
    ) -> nx.Graph:
        """
        Parameters
        ----------
        nodesA : list[int],
            nodes of network A
        nodesB : list[int]
            nodes of network B
        edges : list[tuple[int, int]]
            connecting edges between network A and B
        weights : list[float]
            weights for the edges.

        Returns
        -------
        nx.Graph
            the resulting graph, containing both subgraphs.
        """
        # The initial "weights" are Scores, which need to be translated to weights.
        weights = list(map(lambda x: 1 - x, weights))
        wedges_map = {(e[0], e[1]): w for e, w in zip(edges, weights)}
        wedges = [(e[0], e[1], w) for e, w in zip(edges, weights)]

        # Build Bi-partite graph
        B = nx.Graph()
        B.add_nodes_from(nodesA, bipartite=0)
        B.add_nodes_from(nodesB, bipartite=1)
        B.add_weighted_edges_from(wedges)

        # Matching
        m_edges = [(e[0], e[1], wedges_map[e]) for e in nx.maximal_matching(B)]
        mg = nx.Graph()
        mg.add_nodes_from(nodesA, bipartite=0)
        mg.add_nodes_from(nodesB, bipartite=1)
        mg.add_weighted_edges_from(m_edges)

        return mg
