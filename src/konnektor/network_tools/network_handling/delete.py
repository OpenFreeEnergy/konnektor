from typing import Union

from gufe import LigandNetwork, LigandAtomMapping, Component


def delete_transformation(network: LigandNetwork,
                          transformation: Union[LigandAtomMapping, tuple[
                              Component, Component]]) -> LigandNetwork:
    """
    Remove the desired edge from the network

    Parameters
    ----------
    network: LigandNetwork
    transformation: :Union[LigandAtomMapping, tuple[Component, Component]]

    Returns
    -------
    LigandNetwork
        returns a copy of the ligand network without the removed edge.
    """
    if (isinstance(transformation, LigandAtomMapping)):
        transformation = (transformation.componentA, transformation.componentB)

    f = lambda m: len({m.componentA, m.componentB}.union(transformation)) != 2
    filtered_edges = list(filter(f, network.edges))

    return LigandNetwork(edges=filtered_edges, nodes=network.nodes)


def delete_component(network: LigandNetwork,
                     component: Component) -> LigandNetwork:
    """
    Remove the desired component, which is a node of the graph, and its
    edges from the network.

    Parameters
    ----------
    network: LigandNetwork
    component: Component
        component to be removed from the network, which is a node of the graph.

    Returns
    -------
    LigandNetwork
        returns a copy of the ligand network without the removed component
         and edges containing the component.
    """
    filtered_nodes = list(filter(lambda n: n != component, network.nodes))

    f = lambda m: component not in (m.componentA, m.componentB)
    filtered_edges = list(filter(f, network.edges))

    return LigandNetwork(edges=filtered_edges, nodes=filtered_nodes)
