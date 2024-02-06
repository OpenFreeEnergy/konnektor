# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

from konnektor import network_planners
from konnektor.network_planners import MaximalNetworkPlanner
from .conf import (atom_mapping_basic_test_files, toluene_vs_others,
                   mol_from_smiles, genScorer,
                   GenAtomMapper, BadMapper, ErrorMapper)

@pytest.mark.parametrize('n_process', [1, 2])
@pytest.mark.parametrize('with_progress', [True, False])
@pytest.mark.parametrize('with_scorer', [True, False])
def test_generate_maximal_network(toluene_vs_others, with_progress,
                                  with_scorer, n_process):
    toluene, others = toluene_vs_others

    mapper = GenAtomMapper()

    scorer = genScorer if with_scorer else None

    planner = MaximalNetworkPlanner(
        mapper=mapper, scorer=scorer, progress=with_progress, nprocesses=n_process)
    network = planner.generate_ligand_network(others + [toluene])

    assert len(network.nodes) == len(others) + 1

    edge_count = len(others) * (len(others) + 1) / 2
    assert len(network.edges) == edge_count

    if scorer:
        for edge in network.edges:
            score = edge.annotations['score']
            assert score == 1.0 / len(edge.componentA_to_componentB)
    else:
        for edge in network.edges:
            assert 'score' not in edge.annotations

