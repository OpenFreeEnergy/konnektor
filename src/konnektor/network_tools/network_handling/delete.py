from gufe import Component, LigandAtomMapping, LigandNetwork


def delete_transformation(
    network: LigandNetwork,  # TODO: rename to ligand_network
    transformation: LigandAtomMapping
    | tuple[Component, Component]
    | list[LigandAtomMapping],  # TODO: rename to edges
    must_stay_connected: bool = True,
) -> LigandNetwork:
    """Remove `transformation` from `network`

    Parameters
    ----------
    network : LigandNetwork
    transformation : LigandAtomMapping | tuple[Component, Component] | list[LigandAtomMapping]
    must_stay_connected : bool, optional
        Require that the resulting network remain connected, by default True,

    Returns
    -------
    LigandNetwork
        Copy of `network` without `transformation`.

    Raises
    ------
    RuntimeError
        If the resulting LigandNetwork is not connected.

    """
    if isinstance(transformation, LigandAtomMapping):
        transformations = [(transformation.componentA, transformation.componentB)]
    elif isinstance(transformation, list) and all(
        isinstance(t, LigandAtomMapping) for t in transformation
    ):
        transformations = [(t.componentA, t.componentB) for t in transformation]

    f = lambda m: not any(len({m.componentA, m.componentB}.union(t)) == 2 for t in transformations)
    filtered_edges = list(filter(f, network.edges))

    new_network = LigandNetwork(edges=filtered_edges, nodes=network.nodes)

    if must_stay_connected and not new_network.is_connected():
        raise RuntimeError("Resulting network is not connected anymore!")

    return new_network


def delete_component(
    network: LigandNetwork,  # TODO: rename to ligand_network
    component: Component | list[Component],
    must_stay_connected: bool = True,
) -> LigandNetwork:
    """Remove the `component` and it corresponding edges from the network.

    Parameters
    ----------
    network : LigandNetwork
    component: Component | list[Component]
    must_stay_connected : bool, optional
        Require that the resulting LigandNetwork is still connected, by default True.

    Returns
    -------
    LigandNetwork
        Copy of `network` without `component` or the corresponding edges.

    Raises
    ------
    RuntimeError
        If the resulting LigandNetwork is not connected.
    """
    if isinstance(component, Component):
        components = [component]

    filtered_nodes = list(filter(lambda n: any(n != c for c in components), network.nodes))

    f = lambda m: any(c not in (m.componentA, m.componentB) for c in components)
    filtered_edges = list(filter(f, network.edges))

    new_network = LigandNetwork(edges=filtered_edges, nodes=filtered_nodes)
    if must_stay_connected and not new_network.is_connected():
        raise RuntimeError("Resulting network is not connected anymore!")

    return new_network
