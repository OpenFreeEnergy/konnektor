# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from gufe import LigandNetwork

from konnektor.network_analysis import get_is_connected
from konnektor.network_planners import NNodeEdgesNetworkGenerator
from konnektor.tests.network_planners.conf import (
    atom_mapping_basic_test_files,
    genScorer,
    GenAtomMapper)


def test_nedges_network_mappers(atom_mapping_basic_test_files):
    ligands = [atom_mapping_basic_test_files['toluene'],
               atom_mapping_basic_test_files['2-naftanol'],
               ]

    mapper = GenAtomMapper()
    planner = NNodeEdgesNetworkGenerator(mapper=mapper,
                                         scorer=genScorer,
                                         target_node_connectivity=2)
    network = planner.generate_ligand_network(components=ligands)

    assert isinstance(network, LigandNetwork)
    assert len(network.nodes) == len(ligands)
    assert len(network.edges) <= len(ligands) * 2
    assert get_is_connected(network)
