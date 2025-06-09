# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor


import networkx as nx

from ._abstract_network_algorithm import _AbstractNetworkAlgorithm


class MstNetworkAlgorithm(_AbstractNetworkAlgorithm):
    def generate_network(self, initial_network: nx.Graph, n_edges: int = None) -> nx.Graph:
        # Flip network scores so we can use minimal algorithm
        # TODO: call `update` to make this lighter, or just use max spanning tree.
        g2 = nx.MultiGraph()
        for e1, e2, d in initial_network.graph.edges(data=True):
            g2.add_edge(e1, e2, weight=-d["score"], object=d["object"])

        min_edges = nx.minimum_spanning_edges(g2, weight="weight", keys=True, data=True)

        nodes = initial_network.nodes
        if n_edges is None:
            n_edges = len(nodes) - 1  # max number of MST edges

        selected_edges = [
            (e1, e2, edge_data)
            for i, (e1, e2, _, edge_data) in enumerate(min_edges)
            if (i < n_edges)
        ]

        mg = nx.Graph()
        mg.add_nodes_from(nodes)
        mg.add_edges_from(ebunch_to_add=selected_edges)
        mg.connected = nx.is_connected(mg)

        return mg
