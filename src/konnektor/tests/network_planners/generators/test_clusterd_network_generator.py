# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import numpy as np
from gufe import LigandNetwork
from sklearn.cluster import KMeans
from konnektor.network_analysis import get_is_connected, get_graph_score
from konnektor.network_planners.generators.clustered_network_generator import (
    ClusteredNetworkGenerator,
)
from konnektor.network_tools.clustering.component_diversity_clustering import (
    ComponentsDiversityClusterer,
)
from konnektor.utils.toy_data import build_random_dataset


def test_clustered_network_planner():
    n_compounds = 40
    components, genMapper, genScorer = build_random_dataset(
        n_compounds=n_compounds, rand_seed=42
    )

    from konnektor.network_planners import RadialLigandNetworkPlanner, MstConcatenator

    sub_networker = RadialLigandNetworkPlanner(mapper=genMapper, scorer=genScorer)
    concatenator = MstConcatenator(mapper=genMapper, scorer=genScorer)

    clusterer = ComponentsDiversityClusterer(cluster=KMeans(n_clusters=3))

    planner = ClusteredNetworkGenerator(
        sub_network_planners=sub_networker, concatenator=concatenator, n_processes=1
    )

    ligand_network = planner(components)

    assert isinstance(ligand_network, LigandNetwork)
    assert len(ligand_network.nodes) == n_compounds
    assert len(planner.clusters) == 3
    assert (
        len(ligand_network.edges)
        == 3 * ((n_compounds // 3) - 1) + (3 * concatenator.n_connecting_edges) + 1
    )
    assert get_is_connected(ligand_network)

    np.testing.assert_allclose(get_graph_score(ligand_network), 25.708691, rtol=0.01)
