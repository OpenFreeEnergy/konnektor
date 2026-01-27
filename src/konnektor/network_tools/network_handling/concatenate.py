# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor


from gufe import Component, LigandNetwork

from konnektor.network_planners.concatenators._abstract_network_concatenator import (
    NetworkConcatenator,
)


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
