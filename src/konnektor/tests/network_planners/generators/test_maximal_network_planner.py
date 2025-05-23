# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

from konnektor.network_planners import MaximalNetworkGenerator
from konnektor.tests.network_planners.conf import (
    BadMapper,
    BadMultiMapper,
    ErrorMapper,
    GenAtomMapper,
    SuperBadMapper,
    genScorer,
)


@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
@pytest.mark.parametrize("with_scorer", [True, False])
def test_generate_maximal_network(toluene_vs_others, with_progress, with_scorer, n_process):
    toluene, others = toluene_vs_others

    mapper = GenAtomMapper()

    scorer = genScorer if with_scorer else None

    planner = MaximalNetworkGenerator(
        mappers=mapper, scorer=scorer, progress=with_progress, n_processes=n_process
    )
    network = planner.generate_ligand_network(others + [toluene])

    assert len(network.nodes) == len(others) + 1

    edge_count = len(others) * (len(others) + 1) / 2
    assert len(network.edges) == edge_count

    if scorer:
        for edge in network.edges:
            score = edge.annotations["score"]
            assert score == 1 - 1.0 / len(edge.componentA_to_componentB)
    else:
        for edge in network.edges:
            assert "score" not in edge.annotations


@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
def test_generate_maximal_network_mapper_error(toluene_vs_others, n_process, with_progress):
    """If no scorer, use the first valid mapper.
    # TODO: do we want this behavior, or should we enforce using the first mapper, then error out?
    """
    toluene, others = toluene_vs_others
    components = others + [toluene]

    planner = MaximalNetworkGenerator(
        mappers=[ErrorMapper(), BadMapper()],
        scorer=None,
        progress=with_progress,
        n_processes=n_process,
    )

    network = planner.generate_ligand_network(components)

    # make sure the BadMapper was used ({0:0})
    assert [e.componentA_to_componentB for e in network.edges] == len(network.edges) * [{0: 0}]


@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
def test_generate_maximal_network_no_scorer(toluene_vs_others, n_process, with_progress):
    """If no scorer is provided, the first mapping of the first mapper should be used."""
    toluene, others = toluene_vs_others
    components = others + [toluene]

    planner = MaximalNetworkGenerator(
        mappers=[BadMultiMapper(), SuperBadMapper(), BadMapper()],
        scorer=None,
        progress=with_progress,
        n_processes=n_process,
    )
    # with pytest.warns(match="Only the first mapper provided will be used: <BadMulti"):
    # TODO: warning isn't working with multiprocessing

    network = planner.generate_ligand_network(components)

    # it should use the mapping ({0:2}) of the first mapper (BadMultiMapper)
    assert [e.componentA_to_componentB for e in network.edges] == len(network.edges) * [{0: 2}]
