# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import numpy as np
from gufe import LigandNetwork
from sklearn.cluster import KMeans

from konnektor.network_analysis import get_is_connected, get_network_score
from konnektor.network_planners.generators.clustered_network_generator import (
    StarrySkyNetworkGenerator,
)
from konnektor.network_tools.clustering.component_diversity_clustering import (
    ComponentsDiversityClusterer,
)
from konnektor.utils.toy_data import build_random_dataset


def test_starry_sky_network_planner():
    n_compounds = 40
    components, empty_mapper, random_scorer = build_random_dataset(
        n_compounds=n_compounds, rand_seed=42
    )
    clusterer = ComponentsDiversityClusterer(cluster=KMeans(n_clusters=3))

    planner = StarrySkyNetworkGenerator(
        mappers=empty_mapper, scorer=random_scorer, clusterer=clusterer
    )

    ligand_network = planner(components)
    n_clusters = len(planner.clusters)
    n_connecting_edges = 2
    approx_edges = n_clusters * ((n_compounds // n_clusters) - 1) + (
        n_connecting_edges * n_clusters
    )
    assert isinstance(ligand_network, LigandNetwork)
    assert len(ligand_network.nodes) == n_compounds
    np.testing.assert_allclose(actual=len(ligand_network.edges), desired=approx_edges, rtol=5)
    assert get_is_connected(ligand_network)

    np.testing.assert_allclose(get_network_score(ligand_network), 24.607684, rtol=0.01)
