# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from multiprocessing import Value
import pytest

from konnektor.network_planners import MaximalNetworkGenerator
from konnektor.tests.network_planners.conf import (
    BadMapper,
    ErrorMapper,
    GenAtomMapper,
    MultiAtomMapper,
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

    n_expected_nodes = len(others) + 1
    assert len(network.nodes) == n_expected_nodes

    n_expected_nodes = n_expected_nodes * (n_expected_nodes - 1) / 2
    assert len(network.edges) == n_expected_nodes

    if scorer:
        for edge in network.edges:
            score = edge.annotations["score"]
            assert score == 1.0 / len(edge.componentA_to_componentB)
    else:
        for edge in network.edges:
            assert "score" not in edge.annotations


# TODO: add test for maximal network with missing edges

@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
def test_generate_maximal_network_mapper_error(toluene_vs_others, n_process, with_progress):
    """ """

    toluene, others = toluene_vs_others
    components = others + [toluene]

    # TODO: check that just one error returns an incomplete network
    planner = MaximalNetworkGenerator(
            mappers=[ErrorMapper()],
            scorer=None,
            progress=with_progress,
            n_processes=n_process,
        )
    with pytest.raises(RuntimeError):
        planner.generate_ligand_network(components)



@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
def test_generate_maximal_network_no_scorer_single_mapper(
    toluene_vs_others, n_process, with_progress
):
    """
    """

    toluene, others = toluene_vs_others
    components = others + [toluene]

    planner = MaximalNetworkGenerator(
        mappers=[BadMapper()],
        scorer=None,
        progress=with_progress,
        n_processes=n_process,
    )

    network = planner.generate_ligand_network(components)

    assert [e.componentA_to_componentB for e in network.edges] == len(network.edges) * [{0: 0}]


@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
def test_generate_maximal_network_no_scorer_multiple_mappers(n_process, with_progress):
    """If no scorer but multiple mappers are provided, should throw an error."""

    with pytest.raises(ValueError, match="You must provide a scorer"):
        MaximalNetworkGenerator(
            mappers=[SuperBadMapper(), GenAtomMapper(), BadMapper()],
            scorer=None,
            progress=with_progress,
            n_processes=n_process,
        )

@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
def test_generate_maximal_network_no_scorer_multiple_mappings(
    toluene_vs_others, n_process, with_progress
):
    """
    """

    toluene, others = toluene_vs_others
    components = others + [toluene]
    planner = MaximalNetworkGenerator(
            mappers=[MultiAtomMapper()],
            scorer=None,
            progress=with_progress,
            n_processes=n_process,
        )
    with pytest.raises(ValueError, match="when using mappers that generate multiple mappings"):
        planner.generate_ligand_network(components)
