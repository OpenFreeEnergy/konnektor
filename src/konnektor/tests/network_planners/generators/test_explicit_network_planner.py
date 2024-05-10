# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools

from gufe import LigandNetwork

from konnektor.network_analysis import get_is_connected
from konnektor.network_planners.generators.explicit_network_generator import \
    ExplicitNetworkGenerator
from konnektor.utils.toy_data import build_random_dataset


def test_explicit_network_planner():
    n_compounds = 20
    components, genMapper, genScorer = build_random_dataset(
        n_compounds=n_compounds)

    planner = ExplicitNetworkGenerator(genMapper, genScorer, n_processes=1)

    edges = list(itertools.combinations(components, 2))

    ligand_network = planner(edges)

    assert isinstance(ligand_network, LigandNetwork)
    assert len(ligand_network.nodes) == n_compounds
    assert len(ligand_network.edges) == (n_compounds * (n_compounds - 1)) // 2
    assert get_is_connected(ligand_network)
