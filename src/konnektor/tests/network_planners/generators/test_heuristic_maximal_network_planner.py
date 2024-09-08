# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

from konnektor.network_analysis import get_is_connected
from konnektor.network_planners import HeuristicMaximalNetworkGenerator
from konnektor.utils.toy_data import build_random_dataset


@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
@pytest.mark.parametrize("with_scorer", [True, False])
def test_generate_maximal_network(with_progress, with_scorer, n_process):
    n_compounds = 20
    components, genMapper, genScorer = build_random_dataset(n_compounds=n_compounds)

    planner = HeuristicMaximalNetworkGenerator(
        mappers=genMapper,
        scorer=genScorer,
        n_samples=10,
        progress=with_progress,
        n_processes=n_process,
    )
    network = planner.generate_ligand_network(components)

    assert len(network.nodes) == n_compounds

    edge_count = n_compounds * 10
    assert len(network.edges) <= edge_count
    assert len(network.edges) > n_compounds
    assert get_is_connected(network)
