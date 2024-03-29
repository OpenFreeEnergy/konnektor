import logging
import numpy as np

from sklearn.base import TransformerMixin, ClusterMixin
from sklearn.cluster import KMeans
from sklearn.pipeline import Pipeline
from scikit_mol.fingerprints import RDKitFingerprintTransformer, MorganFingerprintTransformer

from gufe import SmallMoleculeComponent

log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)


class CompoundDiversityClustering():
    def __init__(self, featurize:TransformerMixin = RDKitFingerprintTransformer(),
                cluster:ClusterMixin = KMeans(n_clusters=5, n_init="auto")):
        self.featurize = featurize
        self.cluster = cluster

    def cluster_compounds(self, compounds:list[SmallMoleculeComponent]) -> dict[
        int, SmallMoleculeComponent]:
        # Build Pipeline
        self.pipe = Pipeline([('mol_transformer', self.featurize), ('Cluster', self.cluster)])
        self.pipe.fit([c.to_rdkit() for c in compounds])
        # Retrieve Results
        labels = self.cluster.labels_

        # Compounds label
        cluster_compounds = {}
        for clusterID in np.unique(labels):
            cluster_compounds[clusterID] = [compounds[i] for i,l in enumerate(labels) if(l==clusterID)]

        return cluster_compounds