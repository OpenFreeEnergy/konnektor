# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
from gufe import LigandNetwork

from konnektor.network_planners.concatenators import MstConcatenator
from konnektor.utils.toy_data import EmptyMapper, RandomScorer, build_n_random_mst_network


@pytest.mark.parametrize("n_sub_networks", [2, 3, 4, 6])
def test_mst_concatenation_is_spanning_tree(n_sub_networks):
    n_compounds = 30
    networks = build_n_random_mst_network(
        n_compounds=n_compounds,
        sub_networks=n_sub_networks,
        overlap=0,
        rand_seed=42,
    )
    concatenator = MstConcatenator(
        EmptyMapper(),
        RandomScorer(n=n_compounds),
    )

    connected_network = concatenator.concatenate_networks(ligand_networks=networks)

    # Check the network is connected
    assert connected_network.is_connected()
    # Check we didn't loose any ligands
    assert len(connected_network.nodes) == n_compounds
    # Check that the subnetworks were connected as an MST, meaning k-1 sub-network connections
    n_edges_new = len(connected_network.edges) - sum(len(n.edges) for n in networks)
    assert n_edges_new == n_sub_networks - 1


def test_concatenate_rejects_disconnected_input():
    a, b = build_n_random_mst_network(n_compounds=20, sub_networks=2, overlap=0, rand_seed=42)
    disconnected = LigandNetwork(nodes=a.nodes | b.nodes, edges=a.edges | b.edges)
    concatenator = MstConcatenator(EmptyMapper(), RandomScorer(n=20))
    with pytest.raises(RuntimeError, match="are disconnected"):
        concatenator.concatenate_networks(ligand_networks=[disconnected])


def test_tied_scores_pick_highest_key():
    """On a score tie, the mapping with the max .key is chosen (deterministic)."""
    n = 20
    networkA, networkB = build_n_random_mst_network(
        n_compounds=n, sub_networks=2, overlap=0, rand_seed=42
    )

    class ConstantScorer:
        def __call__(self, mapping):
            return 0.5

    concatenator = MstConcatenator(EmptyMapper(), ConstantScorer())

    # every candidate between the two subnetworks scores 0.5
    candidates = concatenator._score_pair_edges(networkA, networkB)
    expected = max(candidates, key=lambda m: m.key)

    bridges = concatenator._connect_subnetworks_mst([networkA, networkB])
    assert len(bridges) == 1
    assert bridges[0] == expected
