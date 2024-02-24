

from gufe import LigandNetwork, LigandAtomMapping, SmallMoleculeComponent
from ..network_planners.mst_concatenator import MstConcatenate

def merge_two_networks(network1:LigandNetwork,
                    network2:LigandNetwork)->LigandNetwork:
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

    if(len(connecting_nodes)==0):
        raise ValueError("No shared Nodes between Networks")

    merged_nodes = network1.nodes.union(network2.nodes)
    merged_edges = network1.edges.union(network2.edges)

    merged_network = LigandNetwork(edges=merged_edges, nodes=merged_nodes)

    return merged_network

def concatenate(network1:LigandNetwork,
                    network2:LigandNetwork,
           concatenator:MstConcatenate=MstConcatenate)->LigandNetwork:
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
        returns the concatenated network
    """

    concatenated_network = concatenator.concatenate_networks([network1,
                                                             network2])

    return concatenated_network

def append(network:LigandNetwork,
           compound:SmallMoleculeComponent,
           concatenator:MstConcatenate=MstConcatenate)->LigandNetwork:
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
        returns the new network
    """

    appended_network = concatenator.concatenate_networks([network,
                                                             LigandNetwork(
                                                                 edges=[],
                                                                 nodes=[compound])])

    return appended_network

def delete_edge(network:LigandNetwork,
                edge:tuple[SmallMoleculeComponent, SmallMoleculeComponent])->LigandNetwork:

    if(isinstance(edge, LigandAtomMapping)):
        edge = (edge.componentA, edge.componentB)
    
    f = lambda m: len({m.componentA, m.componentB}.union(edge)) != 2
    filtered_edges = filter(f, network.edges)

    return LigandNetwork(edges=filtered_edges, nodes=network.nodes)