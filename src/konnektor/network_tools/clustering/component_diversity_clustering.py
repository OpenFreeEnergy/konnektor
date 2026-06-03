# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import logging

import numpy as np
from gufe import Component
from scikit_mol.fingerprints import MorganFingerprintTransformer
from sklearn.base import ClusterMixin, TransformerMixin
from sklearn.cluster import KMeans
from sklearn.pipeline import Pipeline

from ._abstract_clusterer import _AbstractClusterer

log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)


class ComponentsDiversityClusterer(_AbstractClusterer):
    def __init__(
        self,
        featurize: TransformerMixin = MorganFingerprintTransformer(),  # TODO: move this instantiation into init?
        cluster: ClusterMixin = KMeans(n_clusters=5, n_init="auto"),
        n_processes: int = 1,
    ):
        """This class can be used to separate components by different features, like charge or morgan fingerprints.

        Parameters
        ----------
        featurize : TransformerMixin, optional
            A scikit-learn and scikit-mol compatible featurizer, takes a rdkit mol and transforms to an np.array[number].
            By default MorganFingerprintTransformer().
        cluster : ClusterMixin, optional
            Clustering algorithm compatible with scikit-learn, by default KMeans(n_clusters=5, n_init="auto")
        n_processes : int, optional
            Number of processes that can be used for the network generation, by default 1.
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
        """Featurize and cluster `components` according to their features.

        Parameters
        ----------
        components : list[Component]

        Returns
        -------
        dict[int, list[Component]]
            Clustered compounds, represented as {`clusterid`: [Component]}
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
