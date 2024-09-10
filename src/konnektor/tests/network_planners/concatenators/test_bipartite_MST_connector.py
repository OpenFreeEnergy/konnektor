# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from gufe import LigandNetwork

from konnektor.network_analysis import get_is_connected
from konnektor.network_planners.concatenators import MstConcatenator
from konnektor.tests.network_planners.conf import (
    GenAtomMapper,
    genScorer,
    ligand_network_ab,
    atom_mapping_basic_test_files,
)


# more test here also for the params
def test_mst_network_concatenation(ligand_network_ab):
    concatenator = MstConcatenator(mappers=GenAtomMapper(), scorer=genScorer)

    ln_a, ln_b = ligand_network_ab
    nA = len(ln_a.nodes)
    nB = len(ln_b.nodes)
    eA = len(ln_a.edges)
    eB = len(ln_b.edges)

    cn = concatenator.concatenate_networks([ln_a, ln_b])

    assert isinstance(cn, LigandNetwork)
    assert len(cn.nodes) == nA + nB
    assert len(cn.edges) == eA + eB + concatenator.n_connecting_edges
    assert get_is_connected(cn)
