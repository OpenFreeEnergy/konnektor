# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools

import pytest
from gufe import LigandNetwork

from konnektor.network_analysis import get_is_connected
from konnektor.network_planners.generators.explicit_network_generator import (
    ExplicitNetworkGenerator,
)
from konnektor.utils.toy_data import build_random_dataset


def test_explicit_network_planner():
    n_compounds = 20
    components, emptyMapper, genScorer = build_random_dataset(n_compounds=n_compounds)

    planner = ExplicitNetworkGenerator(emptyMapper, genScorer, n_processes=1)

    edges = list(itertools.combinations(components, 2))

    ligand_network = planner(edges)

    assert isinstance(ligand_network, LigandNetwork)
    assert ligand_network.nodes == frozenset(components)

    result_edges = [(e.componentA, e.componentB) for e in ligand_network.edges]
    assert frozenset(result_edges) == frozenset(edges)

    assert get_is_connected(ligand_network)


def test_explicit_network_planner_from_indices():
    n_compounds = 6
    components, emptyMapper, genScorer = build_random_dataset(n_compounds=n_compounds)
    edges = [(0, 2), (1, 2), (2, 3), (3, 4), (2, 5)]

    planner = ExplicitNetworkGenerator(emptyMapper, genScorer, n_processes=1)

    ligand_network = planner.generate_network_from_indices(components=components, indices=edges)

    assert isinstance(ligand_network, LigandNetwork)
    assert ligand_network.nodes == frozenset(components)

    # names are equal to their indices in this dataset,
    # so we can verify we pulled the correct indices this way:
    result_edges = [(int(e.componentA.name), int(e.componentB.name)) for e in ligand_network.edges]
    assert frozenset(result_edges) == frozenset(edges)

    assert get_is_connected(ligand_network)


def test_explicit_network_planner_from_indices_bad_index():
    n_compounds = 4
    components, emptyMapper, genScorer = build_random_dataset(n_compounds=n_compounds)
    planner = ExplicitNetworkGenerator(emptyMapper, genScorer, n_processes=1)
    with pytest.raises(IndexError, match=r"Invalid ligand index. Requested \(0, 4\)"):
        planner.generate_network_from_indices(components=components, indices=[(0, 4)])


def test_explicit_network_planner_from_indices_disconnected():
    n_compounds = 20
    components, emptyMapper, genScorer = build_random_dataset(n_compounds=n_compounds)
    edges = [(1, 2), (2, 3), (3, 4), (2, 5), (2, 6)]

    planner = ExplicitNetworkGenerator(emptyMapper, genScorer, n_processes=1)

    with pytest.warns(match="Generated network is not connected"):
        ligand_network = planner.generate_network_from_indices(components=components, indices=edges)

    assert isinstance(ligand_network, LigandNetwork)
    assert ligand_network.nodes == frozenset(components)

    # names are equal to their indices in this dataset,
    # so we can verify we pulled the correct indices this way:
    result_edges = [(int(e.componentA.name), int(e.componentB.name)) for e in ligand_network.edges]
    assert frozenset(result_edges) == frozenset(edges)

    assert not get_is_connected(ligand_network)


def test_explicit_network_planner_from_names():
    n_compounds = 6
    components, emptyMapper, genScorer = build_random_dataset(n_compounds=n_compounds)
    edges = [("0", "1"), ("1", "2"), ("2", "3"), ("3", "4"), ("3", "5")]

    planner = ExplicitNetworkGenerator(emptyMapper, genScorer, n_processes=1)

    ligand_network = planner.generate_network_from_names(components=components, names=edges)

    assert isinstance(ligand_network, LigandNetwork)
    assert ligand_network.nodes == frozenset(components)

    result_edges = [(e.componentA.name, e.componentB.name) for e in ligand_network.edges]
    assert frozenset(result_edges) == frozenset(edges)

    assert get_is_connected(ligand_network)


def test_explicit_network_planner_from_names_bad_name():
    n_compounds = 4
    components, emptyMapper, genScorer = build_random_dataset(n_compounds=n_compounds)
    planner = ExplicitNetworkGenerator(emptyMapper, genScorer, n_processes=1)
    with pytest.raises(KeyError, match=r"Invalid name\(s\) requested \['4'\]"):
        planner.generate_network_from_names(components=components, names=[("0", "4")])


def test_explicit_network_planner_from_names_duplicate_name():
    n_compounds = 4
    components, emptyMapper, genScorer = build_random_dataset(n_compounds=n_compounds)
    # make a duplicate
    components.append(components[0])

    planner = ExplicitNetworkGenerator(emptyMapper, genScorer, n_processes=1)
    with pytest.raises(ValueError, match=r"Duplicate names: \['0'\]"):
        planner.generate_network_from_names(components=components, names=[("0", "1")])


def test_explicit_network_planner_from_names_disconnected():
    n_compounds = 20
    components, emptyMapper, genScorer = build_random_dataset(n_compounds=n_compounds)
    edges = [("0", "1"), ("1", "2"), ("2", "3"), ("1", "5"), ("1", "6")]

    planner = ExplicitNetworkGenerator(emptyMapper, genScorer, n_processes=1)

    with pytest.warns(match="Generated network is not connected"):
        ligand_network = planner.generate_network_from_names(components=components, names=edges)

    assert isinstance(ligand_network, LigandNetwork)
    assert ligand_network.nodes == frozenset(components)

    result_edges = [(e.componentA.name, e.componentB.name) for e in ligand_network.edges]
    assert frozenset(result_edges) == frozenset(edges)

    assert not get_is_connected(ligand_network)
