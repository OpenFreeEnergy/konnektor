# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import numpy as np
from gufe import LigandNetwork

from konnektor.network_analysis import get_is_connected, get_network_score
from konnektor.network_planners import NNodeEdgesNetworkGenerator
from konnektor.tests.network_planners.conf import (
    GenAtomMapper,
    atom_mapping_basic_test_files,
    genScorer,
)


def test_nedges_network_mappers(atom_mapping_basic_test_files):
    ligands = [
        atom_mapping_basic_test_files["toluene"],
        atom_mapping_basic_test_files["2-naftanol"],
    ]

    mapper = GenAtomMapper()
    planner = NNodeEdgesNetworkGenerator(
        mappers=mapper, scorer=genScorer, target_component_connectivity=2
    )
    network = planner.generate_ligand_network(components=ligands)

    assert isinstance(network, LigandNetwork)
    assert len(network.nodes) == len(ligands)
    assert len(network.edges) <= len(ligands) * 2
    assert get_is_connected(network)

    np.testing.assert_allclose(get_network_score(network), 0.066667, rtol=0.01)
