# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import numpy as np
import networkx as nx
import pytest

import gufe
from gufe import LigandNetwork

from konnektor.network_planners import RedundantMinimalSpanningTreeNetworkGenerator
from konnektor.tests.network_planners.conf import (
    toluene_vs_others,
    atom_mapping_basic_test_files,
    mol_from_smiles,
    genScorer,
    GenAtomMapper,
    ErrorMapper,
)
from konnektor.network_analysis import get_graph_score


def test_minimal_spanning_network_mappers(atom_mapping_basic_test_files):
    ligands = [
        atom_mapping_basic_test_files["toluene"],
        atom_mapping_basic_test_files["2-naftanol"],
    ]

    mapper = GenAtomMapper()
    planner = RedundantMinimalSpanningTreeNetworkGenerator(
        mapper=mapper, scorer=genScorer, n_redundancy=1
    )
    network = planner.generate_ligand_network(components=ligands)

    assert isinstance(network, LigandNetwork)
    assert list(network.edges)
    np.testing.assert_allclose(get_graph_score(network), 0.066667, rtol=0.01)


@pytest.fixture(scope="session")
def rminimal_spanning_network_redundancy(toluene_vs_others):
    toluene, others = toluene_vs_others
    mapper = GenAtomMapper()
    nred = 3
    planner = RedundantMinimalSpanningTreeNetworkGenerator(
        mapper=mapper, scorer=genScorer, n_redundancy=nred
    )
    network = planner.generate_ligand_network(components=others + [toluene])

    return network, nred


def test_minimal_spanning_network(
    rminimal_spanning_network_redundancy, toluene_vs_others
):
    tol, others = toluene_vs_others
    minimal_spanning_network = rminimal_spanning_network_redundancy[0]
    assert len(minimal_spanning_network.nodes) == len(others) + 1
    for edge in minimal_spanning_network.edges:
        assert edge.componentA_to_componentB != {0: 0}  # lomap should find something


def test_minimal_spanning_network_connectedness(rminimal_spanning_network_redundancy):
    minimal_spanning_network = rminimal_spanning_network_redundancy[0]
    found_pairs = set()
    for edge in minimal_spanning_network.edges:
        pair = frozenset([edge.componentA, edge.componentB])
        assert pair not in found_pairs
        found_pairs.add(pair)

    assert nx.is_connected(nx.MultiGraph(minimal_spanning_network.graph))


@pytest.mark.skip(reason="stochastically wrong and right")
def test_minimal_spanning_network_regression(rminimal_spanning_network_redundancy):
    # this is stochastically failing and working.!
    nred = rminimal_spanning_network_redundancy[1]
    minimal_spanning_network = rminimal_spanning_network_redundancy[0]

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

    print(len(ref), nred, (len(ref) * nred) - nred, len(edge_ids))
    assert len(edge_ids) == (len(ref) * nred) - nred

    # assert edge_ids == ref #This should not be tested here! go to MST generator


@pytest.mark.skip
def test_minimal_spanning_network_unreachable(toluene_vs_others):
    toluene, others = toluene_vs_others
    nimrod = gufe.SmallMoleculeComponent(mol_from_smiles("N"))

    mapper = ErrorMapper()

    with pytest.raises(RuntimeError, match="Unable to create edges"):
        planner = RedundantMinimalSpanningTreeNetworkGenerator(
            mapper=mapper, scorer=genScorer
        )
        network = planner.generate_ligand_network(components=others + [toluene, nimrod])
