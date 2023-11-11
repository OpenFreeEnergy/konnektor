# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import networkx as nx
from collections import defaultdict

from ._abstract_network_generator import _AbstractNetworkGenerator


class StarrySkyNetworkGenerator(_AbstractNetworkGenerator):

    def __int__(self, target_node_connectivity: int = 2):
        self.target_node_connectivity = target_node_connectivity

    def generate_network(self, edges, weights) -> nx.Graph:
        # Build WEdges
        wedges = []
        for edge, weight in zip(edges, weights):
            wedges.append([edge[0], edge[1], weight])

        # Select Edges sorted by weight
        node_edges = defaultdict(list)
        for e1, e2, w in sorted(wedges, key=lambda x: x[2]):
            if (e1 in node_edges and len(node_edges[e1]) < self.target_node_connectivity):
                node_edges[e1].append((e1, e2, w))

            if (e2 in node_edges and len(node_edges[e2]) < self.target_node_connectivity):
                node_edges[e2].append((e1, e2, w))

        # Return Final Graph
        self.g = nx.Graph()
        self.g.add_weighted_edges_from(ebunch_to_add=node_edges)

        return self.g
