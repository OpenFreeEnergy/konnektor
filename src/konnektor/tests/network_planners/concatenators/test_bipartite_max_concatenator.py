# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
from gufe import LigandNetwork

from konnektor.network_planners.concatenators.max_concatenator import MaxConcatenator
from konnektor.tests.network_planners.conf import (
    GenAtomMapper,
    length_scorer,
)
from konnektor.utils.toy_data import EmptyMapper, RandomScorer, build_n_random_mst_network


# more test here also for the params
@pytest.mark.parametrize("n_process", [1, 2])
def test_max_network_concatenation(ligand_network_ab, n_process):
    concatenator = MaxConcatenator(
        mappers=GenAtomMapper(), scorer=length_scorer, n_processes=n_process
    )

    ln_a, ln_b = ligand_network_ab
    nA = len(ln_a.nodes)
    nB = len(ln_b.nodes)
    eA = len(ln_a.edges)
    eB = len(ln_b.edges)

    cn = concatenator.concatenate_networks([ln_a, ln_b])

    assert isinstance(cn, LigandNetwork)
    assert len(cn.nodes) == nA + nB
    assert len(cn.edges) == eA + eB + nA * nB


def test_avoid_edges_excludes_candidate():
    """The avoid_edges should never get scored."""
    n_compounds = 20
    networkA, networkB = build_n_random_mst_network(
        n_compounds=n_compounds,
        sub_networks=2,
        overlap=0,
        rand_seed=42,
    )
    networks = [networkA, networkB]

    # First don't exclude any edges
    concatenator = MaxConcatenator(
        EmptyMapper(),
        RandomScorer(n=n_compounds),
    )
    base = concatenator.concatenate_networks(networks)
    original = set().union(*(n.edges for n in networks))
    avoided = [e for e in base.edges if e not in original][:2]
    # Re-run with those mappings excluded.
    result = concatenator.concatenate_networks(networks, avoid_edges=avoided)

    assert avoided not in list(result.edges)