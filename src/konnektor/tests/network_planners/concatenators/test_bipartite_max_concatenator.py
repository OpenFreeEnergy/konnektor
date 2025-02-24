# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
from gufe import LigandNetwork

from konnektor.network_planners.concatenators.max_concatenator import MaxConcatenator
from konnektor.tests.network_planners.conf import (
    GenAtomMapper,
    genScorer,
)


# more test here also for the params
@pytest.mark.parametrize("n_process", [1, 2])
def test_max_network_concatenation(ligand_network_ab, n_process):
    concatenator = MaxConcatenator(mappers=GenAtomMapper(), scorer=genScorer, n_processes=n_process)

    ln_a, ln_b = ligand_network_ab
    nA = len(ln_a.nodes)
    nB = len(ln_b.nodes)
    eA = len(ln_a.edges)
    eB = len(ln_b.edges)

    cn = concatenator.concatenate_networks([ln_a, ln_b])

    assert isinstance(cn, LigandNetwork)
    assert len(cn.nodes) == nA + nB
    assert len(cn.edges) == eA + eB + nA * nB
