# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import networkx as nx

from ._abstract_network_algorithm import _AbstractNetworkAlgorithm


class MstNetworkAlgorithm(_AbstractNetworkAlgorithm):

    def generate_network(self, edges: list[tuple[int, int]],
                         weights: list[float], n_edges:int=None) -> nx.Graph:
        wedges = []
        nodes = []
        for edge, weight in zip(edges, weights):
            wedges.append([edge[0], edge[1], weight])
            nodes.extend(list(edge))

        if (n_edges is None):
            n_edges = len(nodes) - 1  # max number of MST edges

        self.g = nx.Graph()
        self.g.add_weighted_edges_from(ebunch_to_add=wedges)

        # Next analyze that network to create minimal spanning network. Because
        # we carry the original (directed) AtomMapping, we don't lose
        # direction information when converting to an undirected graph.
        min_edges = nx.minimum_spanning_edges(self.g, weight='weight')
        mse = [(e1, e2, edge_data['weight']) for i, (e1, e2, edge_data) in
               enumerate(min_edges) if (i < n_edges)]

        mg = nx.Graph()
        mg.add_nodes_from(nodes)
        mg.add_weighted_edges_from(ebunch_to_add=mse)
        mg.connected = nx.is_connected(mg)
        return mg
