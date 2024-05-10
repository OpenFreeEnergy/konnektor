# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

from gufe import LigandNetwork
from konnektor.network_planners.concatenator.max_concatenator import MaxConcatenate
from konnektor.tests.network_planners.conf import (GenAtomMapper, genScorer,
                                                   atom_mapping_basic_test_files, ligand_network_ab)


#more test here also for the params
def test_max_network_concatenation(ligand_network_ab):
    concatenator = MaxConcatenate(mapper=GenAtomMapper(), scorer=genScorer)

    ln_a, ln_b = ligand_network_ab
    nA = len(ln_a.nodes)
    nB = len(ln_b.nodes)
    eA = len(ln_a.edges)
    eB = len(ln_b.edges)

    cn = concatenator.concatenate_networks([ln_a, ln_b])

    assert isinstance(cn, LigandNetwork)
    assert len(cn.nodes) == nA + nB
    assert len(cn.edges) == eA + eB + nA * nB