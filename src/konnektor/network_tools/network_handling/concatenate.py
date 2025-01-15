# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from typing import Iterable

from gufe import LigandNetwork, Component

from konnektor.network_planners.concatenators._abstract_network_concatenator import (
    NetworkConcatenator,
)


def concatenate_networks(
    networks: Iterable[LigandNetwork], concatenator: NetworkConcatenator
) -> LigandNetwork:
    """
    Concatenate networks, is similar to a union of a set of nodes and edgees,
    if they are all connected via at least one edge. This means, that  at
    least one node needs to be present in network1 and network2. If this is
    not the case, use the network concatenators.

    Parameters
    ----------
    networks: Iterable[LigandNetwork]
        Network 1 for merging
    concatenator: NetworkConcatenator
        A network planner, that concatenates networks.

    Returns
    -------
    LigandNetwork
        returns the concatenated network
    """

    concat_network = concatenator.concatenate_networks(ligand_networks=networks)

    return concat_network


def append_component(
    network: LigandNetwork, component: Component, concatenator: NetworkConcatenator
) -> LigandNetwork:
    """
    Add one node to the network, based on the provided concatenators algorithm.

    Parameters
    ----------
    network: LigandNetwork
        Network 1 for merging
    component: Component
        append node to network
    concatenator: NetworkConcatenator
        A network planner, that concatenates networks.

    Returns
    -------
    LigandNetwork
        returns the new network
    """
    single_node_net = LigandNetwork(edges=[], nodes=[component])
    appended_network = concatenator.concatenate_networks([network, single_node_net])

    return appended_network
