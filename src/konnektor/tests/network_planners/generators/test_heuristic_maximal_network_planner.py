# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

from konnektor.network_analysis import get_is_connected
from konnektor.network_planners import HeuristicMaximalNetworkGenerator
from konnektor.tests.network_planners.conf import (
    BadMapper,
    BadMultiMapper,
    SuperBadMapper,
)
from konnektor.utils.toy_data import build_random_dataset


@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
def test_generate_maximal_network(with_progress, n_process):
    n_compounds = 20
    components, empty_mapper, random_scorer = build_random_dataset(n_compounds=n_compounds)

    planner = HeuristicMaximalNetworkGenerator(
        mappers=empty_mapper,
        scorer=random_scorer,
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


@pytest.mark.parametrize("with_progress", [True, False])
@pytest.mark.parametrize("n_process", [1, 2])
def test_generate_heuristic_maximal_network_no_scorer(with_progress, n_process):
    """If no scorer is provided, the first mapping of the first mapper should be used."""
    # TODO: what if the first mapper fails? should this be the first *valid* mapper (current behavior), or just error out?
    n_compounds = 4
    components, _, _ = build_random_dataset(n_compounds=n_compounds)

    planner = HeuristicMaximalNetworkGenerator(
        mappers=[BadMultiMapper(), SuperBadMapper(), BadMapper()],
        scorer=None,
        n_samples=3,
        progress=with_progress,
        n_processes=n_process,
    )
    # with pytest.warns(match="Only the first mapper provided will be used: <BadMulti"):
    # TODO: warning isn't working with multiprocessing
    network = planner.generate_ligand_network(components)

    assert len(network.nodes) == n_compounds
    edge_count = n_compounds * 3
    assert len(network.edges) <= edge_count
    assert len(network.edges) > n_compounds
    assert get_is_connected(network)

    # it should use the mapping ({0:2}) of the first mapper (BadMultiMapper)
    assert [e.componentA_to_componentB for e in network.edges] == len(network.edges) * [{0: 2}]
