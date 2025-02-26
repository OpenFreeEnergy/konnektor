from gufe import Component

from konnektor.network_tools.clustering.auxilliary_featurizer import ChargeTransformer
from konnektor.network_tools.clustering.charge_clustering import ChargeClusterer
from konnektor.network_tools.clustering.component_diversity_clustering import (
    ComponentsDiversityClusterer,
)
from konnektor.network_tools.clustering.scaffold_clustering import ScaffoldClusterer
from konnektor.utils.toy_data import build_random_dataset


def test_charge_transformer():
    """only a smoke test! checking that the code runs"""
    n_compounds = 20
    compounds, _, _ = build_random_dataset(n_compounds=n_compounds)

    feat = ChargeTransformer()

    features = feat.fit_transform([c.to_rdkit() for c in compounds])

    assert features.shape[0] == n_compounds
    assert features.shape[1] == 1

    print(features.shape)


def test_diversity_clusterer():
    """only a smoke test! checking that the code runs"""

    n_compounds = 20
    compounds, _, _ = build_random_dataset(n_compounds=n_compounds)

    clusterer = ComponentsDiversityClusterer()
    clusters = clusterer.cluster_compounds(compounds)

    assert isinstance(clusters, dict)
    assert all(isinstance(k, int) for k in clusters.keys())
    assert all(isinstance(v, list) for v in clusters.values())
    assert all(isinstance(c, Component) for v in clusters.values() for c in v)


def test_charge_clusterer():
    """only a smoke test! checking that the code runs"""

    n_compounds = 20
    compounds, _, _ = build_random_dataset(n_compounds=n_compounds)

    clusterer = ChargeClusterer()
    clusters = clusterer.cluster_compounds(compounds)

    assert isinstance(clusters, dict)
    assert all(isinstance(k, int) for k in clusters.keys())
    assert all(isinstance(v, list) for v in clusters.values())
    assert all(isinstance(c, Component) for v in clusters.values() for c in v)


def test_scaffold_clusterer():
    """only a smoke test! checking that the code runs"""

    n_compounds = 20
    compounds, _, _ = build_random_dataset(n_compounds=n_compounds)

    clusterer = ScaffoldClusterer()
    clusters = clusterer.cluster_compounds(compounds)
    print(clusters)
    assert isinstance(clusters, dict)
    assert all(isinstance(k, int) for k in clusters.keys())
    assert all(isinstance(v, list) for v in clusters.values())
    assert all(isinstance(c, Component) for v in clusters.values() for c in v)
