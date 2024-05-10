
from gufe import Component

from konnektor.utils.toy_data import build_random_dataset
from konnektor.network_tools.clustering.charge_clustering import ChargeClusterer
from konnektor.network_tools.clustering.cluster_components import ComponentsDiversityClustering
from konnektor.network_tools.clustering.auxilliary_featurizer import ChargeTransformer


def test_charge_transformer():
    """ only a smoke test! checking that the code runs"""
    n_compounds = 20
    compounds, _, _ = build_random_dataset(n_compounds=n_compounds)

    feat = ChargeTransformer()

    features= feat.fit_transform([c.to_rdkit() for c in compounds])

    assert features.shape[0] == n_compounds
    assert features.shape[1] == 1

    print(features.shape)

def test_charge_clusterer():
    """ only a smoke test! checking that the code runs"""

    n_compounds = 20
    compounds, _, _ = build_random_dataset(n_compounds=n_compounds)

    clusterer = ChargeClusterer()
    clusters = clusterer.cluster_compounds(compounds)

    print(clusters)
    assert isinstance(clusters, dict)
    assert len(clusters) == 1
    assert all(isinstance(k, int) for k in clusters.keys())
    assert all(isinstance(v, list) for v in clusters.values())
    assert all(isinstance(c, Component) for v in clusters.values() for c in v)

def test_diversity_clusterer():
    """ only a smoke test! checking that the code runs"""

    n_compounds = 20
    compounds, _, _ = build_random_dataset(n_compounds=n_compounds)

    clusterer = ComponentsDiversityClustering()
    clusters = clusterer.cluster_compounds(compounds)

    assert isinstance(clusters, dict)
    assert all(isinstance(k, int) for k in clusters.keys())
    assert all(isinstance(v, list) for v in clusters.values())
    assert all(isinstance(c, Component) for v in clusters.values() for c in v)
