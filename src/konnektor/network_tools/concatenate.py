from typing import Iterable
from gufe import LigandNetwork, SmallMoleculeComponent

from ..network_planners.concatenator.mst_concatenator import MstConcatenate


def concatenate_networks(networks:Iterable[LigandNetwork],
                         concatenator)->LigandNetwork:
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

    concat_network = concatenator.concatenate_networks(ligand_networks=networks)

    return concat_network


def append_node(network: LigandNetwork,
                compound: SmallMoleculeComponent,
                concatenator: MstConcatenate)->LigandNetwork:
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
                                                              nodes=
                                                                  [compound])])

    return appended_network
