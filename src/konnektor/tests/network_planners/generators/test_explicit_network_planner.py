# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools

from gufe import LigandNetwork

from konnektor.network_analysis import get_is_connected
from konnektor.network_planners.generators.explicit_network_generator import (
    ExplicitNetworkGenerator,
)
from konnektor.utils.toy_data import build_random_dataset


def test_explicit_network_planner():
    n_compounds = 20
    components, genMapper, genScorer = build_random_dataset(n_compounds=n_compounds)

    planner = ExplicitNetworkGenerator(genMapper, genScorer, n_processes=1)

    edges = list(itertools.combinations(components, 2))

    ligand_network = planner(edges)

    assert isinstance(ligand_network, LigandNetwork)
    assert len(ligand_network.nodes) == n_compounds
    assert len(ligand_network.edges) == (n_compounds * (n_compounds - 1)) // 2
    assert get_is_connected(ligand_network)


def test_explicit_network_planner_with_indices():
    n_compounds = 20
    components, genMapper, genScorer = build_random_dataset(n_compounds=n_compounds)
    indices = [(1, 2), (2, 3), (3, 4), (2, 5), (2, 6)]
    unique_indices = set([i for e in indices for i in e])
    planner = ExplicitNetworkGenerator(genMapper, genScorer, n_processes=1)

    ligand_network = planner.generate_network_from_indices(components=components, indices=indices)

    assert isinstance(ligand_network, LigandNetwork)
    assert len(ligand_network.nodes) == len(unique_indices)
    assert len(ligand_network.edges) == len(indices)
    assert get_is_connected(ligand_network)


def test_explicit_network_planner_with_names():
    n_compounds = 20
    components, genMapper, genScorer = build_random_dataset(n_compounds=n_compounds)
    print(components[0].name)
    names = [("0", "1"), ("1", "2"), ("2", "3"), ("1", "5"), ("1", "6")]
    unique_names = set([i for e in names for i in e])

    planner = ExplicitNetworkGenerator(genMapper, genScorer, n_processes=1)

    ligand_network = planner.generate_network_from_names(components=components, names=names)

    assert isinstance(ligand_network, LigandNetwork)
    assert len(ligand_network.nodes) == len(unique_names)
    assert len(ligand_network.edges) == len(names)
    assert get_is_connected(ligand_network)
