from typing import Iterable

from sklearn.base import TransformerMixin, ClusterMixin
from sklearn.cluster import KMeans
from scikit_mol.fingerprints import RDKitFingerprintTransformer, MorganFingerprintTransformer

from gufe import LigandNetwork, LigandAtomMapping, SmallMoleculeComponent

from .cluster_components import ComponentsDiversityClustering


def delete_transformation(network :LigandNetwork,
                          edge :tuple[SmallMoleculeComponent, SmallMoleculeComponent])\
        ->LigandNetwork:
    if(isinstance(edge, LigandAtomMapping)):
        edge = (edge.componentA, edge.componentB)

    f = lambda m: len({m.componentA, m.componentB}.union(edge)) != 2
    filtered_edges = filter(f, network.edges)

    return LigandNetwork(edges=filtered_edges, nodes=network.nodes)


def cluster_compound(compounds: Iterable[SmallMoleculeComponent],
                     featurize:TransformerMixin = RDKitFingerprintTransformer(),
                    cluster:ClusterMixin = KMeans(n_clusters=5, n_init="auto")
                     )->dict[int,
list[SmallMoleculeComponent]]:
    
    clusterer = ComponentsDiversityClustering(featurize = featurize, cluster=cluster)
    return clusterer.cluster_compounds(compounds)


def cyclize_around_compound(network:LigandNetwork,
                        node:SmallMoleculeComponent)->LigandNetwork:
    #TODO: Implement this
    raise NotImplementedError()