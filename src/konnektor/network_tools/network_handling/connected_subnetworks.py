# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import networkx as nx
from gufe import LigandNetwork


def connected_subnetworks(network: LigandNetwork) -> list[LigandNetwork]:
    """Split a LigandNetwork into its connected subnetworks.

    Parameters
    ----------
    network : LigandNetwork
        The (possibly disconnected) network to split.

    Returns
    -------
    list[LigandNetwork]
        The connected subnetworks of `network`.
    """
    subnetwork_node_sets = list(nx.weakly_connected_components(network.graph))

    if len(subnetwork_node_sets) == 1:
        return [network]

    subnetworks = []
    for subnetwork_nodes in subnetwork_node_sets:
        subnetwork_edges = [
            edge
            for edge in network.edges
            if edge.componentA in component and edge.componentB in subnetwork_nodes
        ]
        subnetworks.append(LigandNetwork(nodes=subnetwork_nodes, edges=subnetwork_edges))

    return subnetworks
