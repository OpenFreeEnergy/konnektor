# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import networkx as nx
from gufe import LigandNetwork


def decompose_network(network: LigandNetwork) -> list[LigandNetwork]:
    """Split a disconnected LigandNetwork into its connected sub-networks.

    Parameters
    ----------
    network : LigandNetwork
        The (possibly disconnected) network to split.

    Returns
    -------
    list[LigandNetwork]
        The connected sub-networks of `network`.
    """
    graph = network.graph.to_undirected()

    sub_networks = []
    for component in nx.connected_components(graph):
        sub_edges = [
            edge
            for edge in network.edges
            if edge.componentA in component and edge.componentB in component
        ]
        sub_networks.append(LigandNetwork(nodes=component, edges=sub_edges))

    return sub_networks
