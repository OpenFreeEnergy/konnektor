# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor


import networkx as nx

from ._abstract_network_algorithm import _AbstractNetworkAlgorithm


class MstNetworkAlgorithm(_AbstractNetworkAlgorithm):
    def generate_network(self, initial_network: nx.Graph, n_edges: int = None) -> nx.Graph:
        # we actually use a maximum spanning tree since higher scores are better
        # TODO: make sure we get the directed edges back out
        min_edges = nx.maximum_spanning_edges(
            initial_network.to_undirected(), weight="score", keys=True, data=True
        )

        nodes = initial_network.nodes(data=True)
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
