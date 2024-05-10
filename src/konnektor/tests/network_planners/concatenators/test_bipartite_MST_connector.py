# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
import numpy as np
import networkx as nx

from gufe import LigandNetwork
from konnektor.network_planners.concatenators import MstConcatenate
from konnektor.tests.network_planners.conf import (GenAtomMapper, genScorer,
                                                   atom_mapping_basic_test_files, ligand_network_ab)


#more test here also for the params
def test_mst_network_concatenation(ligand_network_ab):
    concatenator = MstConcatenate(mapper=GenAtomMapper(), scorer=genScorer)

    ln_a, ln_b = ligand_network_ab

    cn = concatenator.concatenate_networks([ln_a, ln_b])

    assert isinstance(cn, LigandNetwork)
    assert len(cn.nodes) == len(ln_a.nodes) + len(ln_b.nodes)
    assert (len(cn.edges) == len(ln_a.edges) + len(ln_b.edges) +
            concatenator.n_connecting_edges)
