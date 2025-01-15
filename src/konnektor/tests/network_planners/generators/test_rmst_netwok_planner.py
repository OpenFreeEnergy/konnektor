# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import gufe
import networkx as nx
import numpy as np
import pytest
from gufe import LigandNetwork

from konnektor.network_analysis import get_network_score
from konnektor.network_planners import RedundantMinimalSpanningTreeNetworkGenerator
from konnektor.tests.network_planners.conf import (
    ErrorMapper,
    GenAtomMapper,
    atom_mapping_basic_test_files,
    genScorer,
    mol_from_smiles,
    toluene_vs_others,
)


def test_rminimal_spanning_network_mappers(atom_mapping_basic_test_files):
    ligands = [
        atom_mapping_basic_test_files["toluene"],
        atom_mapping_basic_test_files["2-naftanol"],
    ]

    mapper = GenAtomMapper()
    planner = RedundantMinimalSpanningTreeNetworkGenerator(
        mappers=mapper, scorer=genScorer, n_redundancy=1
    )
    network = planner.generate_ligand_network(components=ligands)

    assert isinstance(network, LigandNetwork)
    assert list(network.edges)
    np.testing.assert_allclose(get_network_score(network), 0.066667, rtol=0.01)


@pytest.fixture(scope="session")
def rminimal_spanning_network_redundancy(toluene_vs_others):
    toluene, others = toluene_vs_others
    mapper = GenAtomMapper()
    nred = 3
    planner = RedundantMinimalSpanningTreeNetworkGenerator(
        mappers=mapper, scorer=genScorer, n_redundancy=nred
    )
    network = planner.generate_ligand_network(components=others + [toluene])

    return network, nred


def test_rminimal_spanning_network(rminimal_spanning_network_redundancy, toluene_vs_others):
    tol, others = toluene_vs_others
    minimal_spanning_network = rminimal_spanning_network_redundancy[0]
    assert len(minimal_spanning_network.nodes) == len(others) + 1
    for edge in minimal_spanning_network.edges:
        assert edge.componentA_to_componentB != {0: 0}  # lomap should find something


def test_rminimal_spanning_network_connectedness(rminimal_spanning_network_redundancy):
    minimal_spanning_network = rminimal_spanning_network_redundancy[0]
    found_pairs = set()
    for edge in minimal_spanning_network.edges:
        pair = frozenset([edge.componentA, edge.componentB])
        assert pair not in found_pairs
        found_pairs.add(pair)

    assert nx.is_connected(nx.MultiGraph(minimal_spanning_network.graph))


def test_minimal_rmst_network_noedger(toluene_vs_others):
    toluene, others = toluene_vs_others
    nimrod = gufe.SmallMoleculeComponent(mol_from_smiles("N"))

    mapper = ErrorMapper()

    with pytest.raises(RuntimeError, match="Could not generate any mapping"):
        planner = RedundantMinimalSpanningTreeNetworkGenerator(mappers=mapper, scorer=genScorer)
        network = planner.generate_ligand_network(components=others + [toluene, nimrod])
