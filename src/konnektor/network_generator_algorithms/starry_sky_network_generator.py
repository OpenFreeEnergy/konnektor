# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import networkx as nx
from collections import defaultdict

from ._abstract_network_generator import _AbstractNetworkGenerator


class StarrySkyNetworkGenerator(_AbstractNetworkGenerator):

    def __int__(self, target_node_connectivity: int = 2):
        self.target_node_connectivity = target_node_connectivity

    def generate_network(self, edges, weights) -> nx.Graph:
        w_edges = []
        nodes = []
        for e, w in zip(edges, weights):
            w_edges.append((e[0], e[1], w))
            # w_edges.append((e[1], e[0], w))
            nodes.extend(e)

        node_edges = defaultdict(list)
        node_con = {n: 0 for n in nodes}  # early termination crit
        for e1, e2, w in sorted(w_edges, key=lambda x: x[2]):
            print(e1, e2, w)
            if (len(node_edges[e1]) < self.target_node_connectivity):
                node_edges[e1].append((e1, e2, w))
                node_con[e1] += 1

            if (len(node_edges[e2]) < self.target_node_connectivity):
                node_edges[e2].append((e1, e2, w))
                node_con[e2] += 1

            if (all(v == self.target_node_connectivity for v in node_con.values())):
                break

        edges = list(set([e for el in node_edges.values() for e in el]))

        g = nx.Graph()
        g.add_weighted_edges_from(ebunch_to_add=edges)

        return self.g
