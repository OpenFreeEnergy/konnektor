from typing import Union

from gufe import LigandNetwork, LigandAtomMapping, Component


def delete_transformation(
    network: LigandNetwork,
    transformation: Union[
        LigandAtomMapping, tuple[Component, Component], list[LigandAtomMapping]
    ],
) -> LigandNetwork:
    """
    Remove the desired edge from the network

    Parameters
    ----------
    network: LigandNetwork
    transformation: :Union[LigandAtomMapping, tuple[Component, Component], list[LigandAtomMapping]]

    Returns
    -------
    LigandNetwork
        returns a copy of the ligand network without the removed edge.
    """
    if isinstance(transformation, LigandAtomMapping):
        transformations = [(transformation.componentA, transformation.componentB)]
    elif isinstance(transformation, list) and all(
        isinstance(t, LigandAtomMapping) for t in transformation
    ):
        transformations = [(t.componentA, t.componentB) for t in transformation]

    f = lambda m: not any(
        len({m.componentA, m.componentB}.union(t)) == 2 for t in transformations
    )
    filtered_edges = list(filter(f, network.edges))

    return LigandNetwork(edges=filtered_edges, nodes=network.nodes)


def delete_component(
    network: LigandNetwork, component: Union[Component, list[Component]]
) -> LigandNetwork:
    """
    Remove the desired component, which is a node of the graph, and its
    edges from the network.

    Parameters
    ----------
    network: LigandNetwork
    component:Union[Component, list[Component]]
        component to be removed from the network, which is a node of the graph.

    Returns
    -------
    LigandNetwork
        returns a copy of the ligand network without the removed component
        and edges containing the component.
    """
    if isinstance(component, Component):
        components = [component]

    filtered_nodes = list(
        filter(lambda n: any(n != c for c in components), network.nodes)
    )

    f = lambda m: any(c not in (m.componentA, m.componentB) for c in components)
    filtered_edges = list(filter(f, network.edges))

    return LigandNetwork(edges=filtered_edges, nodes=filtered_nodes)
