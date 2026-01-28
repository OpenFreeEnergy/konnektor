# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import numpy as np
from gufe import LigandNetwork
from sklearn.cluster import KMeans

from konnektor.network_analysis import get_is_connected, get_network_score
from konnektor.network_planners.generators.clustered_network_generator import (
    ClusteredNetworkGenerator,
)
from konnektor.network_tools.clustering.component_diversity_clustering import (
    ComponentsDiversityClusterer,
)
from konnektor.utils.toy_data import build_random_dataset


def test_clustered_network_planner():
    n_compounds = 40

    # TODO: this test is flaky, depending on the random seed
    components, empty_mapper, smiles_length_scorer = build_random_dataset(
        n_compounds=n_compounds, rand_seed=42
    )

    from konnektor.network_planners import MstConcatenator, RadialNetworkGenerator

    sub_networker = RadialNetworkGenerator(mappers=empty_mapper, scorer=smiles_length_scorer)
    concatenator = MstConcatenator(mappers=empty_mapper, scorer=smiles_length_scorer)
    clusterer = ComponentsDiversityClusterer(cluster=KMeans(n_clusters=3))

    planner = ClusteredNetworkGenerator(
        sub_network_planners=sub_networker,
        concatenator=concatenator,
        n_processes=1,
        clusterer=clusterer,
    )

    ligand_network = planner(components)
    assert isinstance(ligand_network, LigandNetwork)
    assert len(ligand_network.nodes) == n_compounds
    assert len(planner.clusters) == 3
    expected_number_of_edges = 3 * ((n_compounds // 3) - 1) + (3 * concatenator.n_connecting_edges)

    assert len(ligand_network.edges) == expected_number_of_edges
    assert get_is_connected(ligand_network)

    np.testing.assert_allclose(get_network_score(ligand_network), 25.708691, rtol=0.05)
