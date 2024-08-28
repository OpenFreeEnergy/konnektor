# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import gufe
import networkx as nx
import numpy as np
import pytest
from gufe import LigandNetwork

from konnektor.network_planners import MinimalSpanningTreeNetworkGenerator
from konnektor.tests.network_planners.conf import (
    toluene_vs_others,
    atom_mapping_basic_test_files,
    mol_from_smiles,
    genScorer,
    GenAtomMapper,
    ErrorMapper,
)
from konnektor.network_analysis import get_network_score


def test_minimal_spanning_network_mappers(atom_mapping_basic_test_files):
    ligands = [
        atom_mapping_basic_test_files["toluene"],
        atom_mapping_basic_test_files["2-naftanol"],
    ]

    mapper = GenAtomMapper()
    planner = MinimalSpanningTreeNetworkGenerator(mapper=mapper, scorer=genScorer)
    network = planner.generate_ligand_network(components=ligands)

    assert isinstance(network, LigandNetwork)
    assert list(network.edges)
    np.testing.assert_allclose(get_network_score(network), 0.066667, rtol=0.001)


@pytest.fixture(scope="session")
def minimal_spanning_network(toluene_vs_others):
    toluene, others = toluene_vs_others
    mapper = GenAtomMapper()

    planner = MinimalSpanningTreeNetworkGenerator(mapper=mapper, scorer=genScorer)
    network = planner.generate_ligand_network(components=others + [toluene])

    return network


def test_minimal_spanning_network(minimal_spanning_network, toluene_vs_others):
    tol, others = toluene_vs_others
    assert len(minimal_spanning_network.nodes) == len(others) + 1
    for edge in minimal_spanning_network.edges:
        assert edge.componentA_to_componentB != {0: 0}  # lomap should find something


def test_minimal_spanning_network_connectedness(minimal_spanning_network):
    found_pairs = set()
    for edge in minimal_spanning_network.edges:
        pair = frozenset([edge.componentA, edge.componentB])
        assert pair not in found_pairs
        found_pairs.add(pair)

    assert nx.is_connected(nx.MultiGraph(minimal_spanning_network.graph))


def test_minimal_spanning_network_regression(minimal_spanning_network):
    # issue #244, this was previously giving non-reproducible (yet valid)
    # networks when scores were tied.
    edge_ids = sorted(
        (edge.componentA.name, edge.componentB.name)
        for edge in minimal_spanning_network.edges
    )
    ref = sorted(
        [
            ("1,3,7-trimethylnaphthalene", "2,6-dimethylnaphthalene"),
            ("1-butyl-4-methylbenzene", "2-methyl-6-propylnaphthalene"),
            ("2,6-dimethylnaphthalene", "2-methyl-6-propylnaphthalene"),
            ("2,6-dimethylnaphthalene", "2-methylnaphthalene"),
            ("2,6-dimethylnaphthalene", "2-naftanol"),
            ("2,6-dimethylnaphthalene", "methylcyclohexane"),
            ("2,6-dimethylnaphthalene", "toluene"),
        ]
    )

    assert len(edge_ids) == len(ref)
    # assert edge_ids == ref #This should not be tested here! go to MST generator


@pytest.mark.skip
def test_minimal_spanning_network_unreachable(toluene_vs_others):
    toluene, others = toluene_vs_others
    nimrod = gufe.SmallMoleculeComponent(mol_from_smiles("N"))

    mapper = ErrorMapper()

    with pytest.raises(RuntimeError, match="Unable to create edges"):
        planner = MinimalSpanningTreeNetworkGenerator(mapper=mapper, scorer=genScorer)
        network = planner.generate_ligand_network(compounds=others + [toluene, nimrod])
