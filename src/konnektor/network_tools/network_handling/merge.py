# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from gufe import LigandNetwork


def merge_two_networks(network1: LigandNetwork,
                       network2: LigandNetwork) -> LigandNetwork:
    """
    Merging networks, is similar to a union of a set of nodes and edgees,
    if they are all connected via at least one edge. This means, that  at
    least one node needs to be present in network1 and network2. If this is
    not the case, use the network concatenators.

    Parameters
    ----------
    network1: LigandNetwork
        Network 1 for merging
    network2
        Network 1 for merging

    Returns
    -------
    LigandNetwork
        returns the merged network
    """
    connecting_nodes = network1.nodes.intersection(network2.nodes)

    if (len(connecting_nodes) == 0):
        raise ValueError("No shared Nodes between Networks")

    merged_nodes = network1.nodes.union(network2.nodes)
    merged_edges = network1.edges.union(network2.edges)
    merged_network = LigandNetwork(edges=merged_edges, nodes=merged_nodes)

    return merged_network


def merge_networks(networks: list[LigandNetwork]) -> LigandNetwork:
    """
    Merging networks, is similar to a union of a set of nodes and edgees,
    if they are all connected via at least one edge. This means, that  at
    least one node needs to be present in network1 and network2. If this is
    not the case, use the network concatenators.

    Todo: find a mergable path in multiple sets of networks

    Parameters
    ----------
    networks: list[LigandNetwork]
        Networks for merging

    Returns
    -------
    LigandNetwork
        returns the merged network
    """

    # merge_path
    # Todo: Implement find a mergeable path

    # do merging
    networkA = networks[0]
    merged_nodes = []
    merged_edges = []
    for networkB in networks[1:]:
        connecting_nodes = networkA.nodes.difference(networkB.nodes)

        if len(connecting_nodes) >= len(networkB.nodes):
            raise ValueError("No shared Nodes between Networks")

        merged_nodes.extend(list(networkA.nodes.union(networkB.nodes)))
        merged_edges.extend(list(networkA.edges.union(networkB.edges)))

    merged_network = LigandNetwork(edges=set(merged_edges),
                                   nodes=set(merged_nodes))

    return merged_network
