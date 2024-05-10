from typing import Iterable
from gufe import LigandNetwork, Component

from ..network_planners.concatenator.mst_concatenator import MstConcatenate


def concatenate_networks(networks: Iterable[LigandNetwork],
                         concatenator) -> LigandNetwork:
    """
    Concatenate networks, is similar to a union of a set of nodes and edgees,
    if they are all connected via at least one edge. This means, that  at
    least one node needs to be present in network1 and network2. If this is
    not the case, use the network concatenators.

    Parameters
    ----------
    networks: Iterable[LigandNetwork]
        Network 1 for merging

    Returns
    -------
    LigandNetwork
        returns the concatenated network
    """

    concat_network = concatenator.concatenate_networks(ligand_networks=networks)

    return concat_network


def append_node(network: LigandNetwork,
                component: Component,
                concatenator: MstConcatenate) -> LigandNetwork:
    """
    Add one node to the network,

    Parameters
    ----------
    network: LigandNetwork
        Network 1 for merging
    component: Component
        append node to network

    Returns
    -------
    LigandNetwork
        returns the new network
    """

    appended_network = concatenator.concatenate_networks([network,
                                                          LigandNetwork(
                                                              edges=[],
                                                              nodes=
                                                              [component])])

    return appended_network
