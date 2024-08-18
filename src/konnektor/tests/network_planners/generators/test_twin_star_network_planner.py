# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
import numpy as np
from gufe import LigandNetwork

from konnektor.network_planners import TwinStarNetworkGenerator
from konnektor.network_analysis import get_is_connected, get_graph_score
from konnektor.utils.toy_data import build_random_dataset


def test_twin_star_network_planner():
    n_compounds = 40
    components, genMapper, genScorer = build_random_dataset(
        n_compounds=n_compounds, rand_seed=42
    )

    planner = TwinStarNetworkGenerator(mapper=genMapper, scorer=genScorer)

    # Testing
    ligand_network = planner(components)
    n_centers = planner.n_centers
    approx_edges = (len(components) - 1) * n_centers
    assert isinstance(ligand_network, LigandNetwork)
    assert len(ligand_network.nodes) == n_compounds
    np.testing.assert_allclose(
        actual=len(ligand_network.edges), desired=approx_edges, rtol=5
    )
    assert get_is_connected(ligand_network)
    np.testing.assert_allclose(get_graph_score(ligand_network), 39.944662, rtol=0.01)
