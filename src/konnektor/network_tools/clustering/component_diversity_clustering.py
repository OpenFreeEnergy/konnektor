# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import logging

import numpy as np
from gufe import Component
from scikit_mol.fingerprints import MorganFingerprintTransformer
from sklearn.base import TransformerMixin, ClusterMixin
from sklearn.cluster import KMeans
from sklearn.pipeline import Pipeline

from ._abstract_clusterer import _AbstractClusterer

log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)


class ComponentsDiversityClusterer(_AbstractClusterer):
    def __init__(
        self,
        featurize: TransformerMixin = MorganFingerprintTransformer(),
        cluster: ClusterMixin = KMeans(n_clusters=5, n_init="auto"),
        n_processes: int = 1,
    ):
        """
        This class can be use to seperate components by different features, like charge or morgan fingerprints.

        Parameters
        ----------
        featurize: TransformerMixin, optional
            A scikit-learn and scikit-mol compatible featurizer, takes a rdkit mol and transforms to an np.array[number].
            As default the morgan fingerprints are used.
        cluster: ClusterMixin
            a scikit-learn compatible clustering algorithm.
            as default a  KMeans(n_clusters=5, n_init="auto") is used.
        parallel: int, optional
            tries to push the parallelization triggers of featurize and cluster

        """
        self._cluster_centers = None
        self.featurize = featurize
        if hasattr(self.featurize, "parallel") and n_processes > 1:
            self.featurize.parallel = n_processes

        self.cluster = cluster
        if hasattr(self.cluster, "n_jobs") and n_processes > 1:
            self.cluster.n_jobs = n_processes

    @property
    def cluster_centers(self) -> int:
        if self._cluster_centers is None:
            raise ValueError("Cluster centers were not set.")
        else:
            return self._cluster_centers

    def cluster_compounds(self, components: list[Component]) -> dict[int, list[Component]]:
        """
            The method featurizes and clusters the molecules according to the features.


        Parameters
        ----------
        components:list[Component]
            the list of components, that should be seperated into different categories.


        Returns
        -------
        dict[int, list[Component]]
            the index represents the clusterid, the values are lists of Components, corresponding to the clusters.
        """
        # Build Pipeline
        self.pipe = Pipeline([("mol_transformer", self.featurize), ("Cluster", self.cluster)])
        self.pipe.fit([c.to_rdkit() for c in components])

        # Retrieve Results
        labels = self.cluster.labels_

        if hasattr(self.cluster, "cluster_centers_"):
            self._cluster_centers = self.cluster.cluster_centers_

        # Compounds label
        cluster_components = {}
        for clusterID in np.unique(labels):
            cluster_components[int(clusterID)] = [
                components[i] for i, cid in enumerate(labels) if (cid == clusterID)
            ]

        return cluster_components
