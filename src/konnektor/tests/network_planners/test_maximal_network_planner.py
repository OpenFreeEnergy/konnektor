# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

import openfe

from konnektor import network_planners

from .conf import atom_mapping_basic_test_files, toluene_vs_others, mol_from_smiles,BadMapper


@pytest.mark.parametrize('with_progress', [True, False])
@pytest.mark.parametrize('with_scorer', [True, False])
def test_generate_maximal_network(toluene_vs_others, with_progress,
                                  with_scorer, extra_mapper):
    toluene, others = toluene_vs_others

    mapper = openfe.setup.atom_mapping.LomapAtomMapper()

    def scoring_func(mapping):
        return 1.0 / len(mapping.componentA_to_componentB)

    scorer = scoring_func if with_scorer else None

    planner = network_planners.MinimalSpanningTreeLigandNetworkPlanner(mapper=mapper, scorer=scoring_func)
    network = planner.generate_ligand_network(ligands=others + [toluene], progress=with_progress)

    assert len(network.nodes) == len(others) + 1

    if extra_mapper:
        edge_count = len(others) * (len(others) + 1)
    else:
        edge_count = len(others) * (len(others) + 1) / 2

    assert len(network.edges) == edge_count

    if scorer:
        for edge in network.edges:
            score = edge.annotations['score']
            assert score == 1.0 / len(edge.componentA_to_componentB)
    else:
        for edge in network.edges:
            assert 'score' not in edge.annotations

