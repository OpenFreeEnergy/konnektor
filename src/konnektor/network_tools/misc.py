from typing import Iterable, Union

from gufe import Component
from gufe import LigandNetwork, LigandAtomMapping

from sklearn.base import TransformerMixin, ClusterMixin
from sklearn.cluster import KMeans
from scikit_mol.fingerprints import RDKitFingerprintTransformer
from konnektor.network_tools.clustering._abstract_clusterer import _AbstractClusterer


def delete_transformation(network: LigandNetwork,
                          edge: Union[LigandAtomMapping, tuple[Component, Component]]) -> LigandNetwork:
    """
    Remove the desired edge from the network

    Parameters
    ----------
    network: LigandNetwork
    edge: :Union[LigandAtomMapping, tuple[Component, Component]]

    Returns
    -------
    LigandNetwork
        returns a copy of the ligand network without the removed edge.
    """
    if (isinstance(edge, LigandAtomMapping)):
        edge = (edge.componentA, edge.componentB)

    f = lambda m: len({m.componentA, m.componentB}.union(edge)) != 2
    filtered_edges = list(filter(f, network.edges))

    return LigandNetwork(edges=filtered_edges, nodes=network.nodes)


def delete_component(network: LigandNetwork,
                     component: Component) -> LigandNetwork:
    """
    Remove the desired component, which is a node of the graph, and its edges from the network.

    Parameters
    ----------
    network: LigandNetwork
    component: Component
        component to be removed from the network, which is a node of the graph.

    Returns
    -------
    LigandNetwork
        returns a copy of the ligand network without the removed component and edges containing the component.
    """
    f = lambda n: n != component
    filtered_nodes = list(filter(f, network.nodes))

    f = lambda m: component not in (m.componentA, m.componentB)
    filtered_edges = list(filter(f, network.edges))

    return LigandNetwork(edges=filtered_edges, nodes=filtered_nodes)


def cluster_components(compounds: Iterable[Component],
                       clusterer:_AbstractClusterer
                       ) -> dict[int, list[Component]]:
    """
    This is a helper function for using the clustering alsgorithms.
    TODO: Implement this
    Not Implemented!

    Parameters
    ----------
    compounds: Iterable[SmallMoleculeComponent]
    clusterer: _AbstractClusterer
        the function used to featurize the

    Returns
    -------

    """
    raise NotImplementedError()


def cyclize_around_component(network: LigandNetwork,
                             node: Component) -> LigandNetwork:
    """
        TODO: Implement this
        Not Implemented!

    Parameters
    ----------
    network
    node

    Returns
    -------

    """
    # TODO: Implement this
    raise NotImplementedError()
