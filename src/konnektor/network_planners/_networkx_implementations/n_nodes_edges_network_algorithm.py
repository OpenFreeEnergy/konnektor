# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import networkx as nx

from ._abstract_network_algorithm import _AbstractNetworkAlgorithm


# Todo: check that the resulting graph is connnected.


class NNodeEdgesNetworkAlgorithm(_AbstractNetworkAlgorithm):
    def __init__(self, target_node_connectivity: int = 2):
        self.target_node_connectivity = target_node_connectivity

    def generate_network(
        self, edges: list[tuple[int, int]], weights: list[float]
    ) -> nx.Graph:
        w_edges = []
        nodes = []
        # The initial "weights" are Scores, which need to be translated to weights.
        weights = list(map(lambda x: 1 - x, weights))
        for e, w in zip(edges, weights):
            w_edges.append((e[0], e[1], w))
            nodes.extend(e)

        # Build initial MST:
        self.g = nx.Graph()
        self.g.add_weighted_edges_from(ebunch_to_add=w_edges)
        min_edges = nx.minimum_spanning_edges(self.g, weight='weight')
        mse = [(e1, e2, edge_data['weight']) for i, (e1, e2, edge_data) in
               enumerate(min_edges)]

        # Select Additional Edges
        final_edges = mse
        rest_edges = list(filter(lambda e: e not in mse, mse))
        node_con = {n: 1 for n in nodes}  # early termination crit
        for edge in sorted(rest_edges, key=lambda x: x[2]):
            e1 = edge[0]
            e2 = edge[1]
            if node_con[e1] < self.target_node_connectivity:
                if edge not in final_edges:
                    final_edges.append(edge)
                node_con[e1] += 1

            if node_con[e2] < self.target_node_connectivity:
                if edge not in final_edges:
                    final_edges.append(edge)
                node_con[e2] += 1

            if all(v == self.target_node_connectivity for v in node_con.values()):
                break

        self.g = nx.Graph()
        self.g.add_weighted_edges_from(ebunch_to_add=final_edges)

        return self.g
