# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
from gufe import LigandNetwork

from konnektor.network_planners.concatenators import MstConcatenator
from konnektor.utils.toy_data import EmptyMapper, RandomScorer, build_n_random_mst_network


@pytest.mark.parametrize("n_sub_networks", [2, 3, 4, 6])
def test_mst_concatenation_is_spanning_tree(n_sub_networks):
    n_compounds = 30
    n_connecting_edges = 1
    networks = build_n_random_mst_network(
        n_compounds=n_compounds,
        sub_networks=n_sub_networks,
        overlap=0,
        rand_seed=42,
    )
    concatenator = MstConcatenator(
        EmptyMapper(),
        RandomScorer(n=n_compounds),
        n_connecting_edges=n_connecting_edges,
    )

    connected_network = concatenator.concatenate_networks(ligand_networks=networks)

    # Check the network is connected
    assert connected_network.is_connected()
    # Check we didn't loose any ligands
    assert len(connected_network.nodes) == n_compounds
    # Check that the subnetworks were connected as an MST, meaning k-1 sub-network connections
    n_edges_new = len(connected_network.edges) - sum(len(n.edges) for n in networks)
    assert n_edges_new == n_connecting_edges * (n_sub_networks - 1)


def test_mst_network_concatenation_redundancy():
    n_connecting_edges = 3
    n_compounds = 30
    n_sub_networks = 4

    networks = build_n_random_mst_network(
        n_compounds=n_compounds, sub_networks=n_sub_networks, overlap=0, rand_seed=42
    )
    concatenator = MstConcatenator(
        EmptyMapper(), RandomScorer(n=n_compounds), n_connecting_edges=n_connecting_edges
    )
    connected_network = concatenator.concatenate_networks(ligand_networks=networks)

    assert connected_network.is_connected()
    n_edges_new = len(connected_network.edges) - sum(len(n.edges) for n in networks)
    assert n_edges_new == n_connecting_edges * (n_sub_networks - 1)


def test_concatenate_rejects_disconnected_input():
    a, b = build_n_random_mst_network(n_compounds=20, sub_networks=2, overlap=0, rand_seed=42)
    disconnected = LigandNetwork(nodes=a.nodes | b.nodes, edges=a.edges | b.edges)
    concatenator = MstConcatenator(EmptyMapper(), RandomScorer(n=20))
    with pytest.raises(RuntimeError, match="are disconnected"):
        concatenator.concatenate_networks(ligand_networks=[disconnected])


def test_avoid_edges_excludes_candidate():
    """The avoid_edges should never get scored."""
    n_compounds = 20
    networkA, networkB = build_n_random_mst_network(
        n_compounds=n_compounds,
        sub_networks=2,
        overlap=0,
        rand_seed=42,
    )

    # First don't exclude any edges
    concatenator = MstConcatenator(
        EmptyMapper(),
        RandomScorer(n=n_compounds),
    )
    mappings = concatenator._score_pair_edges(networkA, networkB)
    avoided = mappings[0]

    # Re-run with that mapping excluded.
    concatenator = MstConcatenator(
        EmptyMapper(),
        RandomScorer(n=n_compounds),
        avoid_edges=[avoided],
    )
    mappings = concatenator._score_pair_edges(networkA, networkB)

    resulting_pairs = {frozenset((mapping.componentA, mapping.componentB)) for mapping in mappings}
    avoided_pair = frozenset((avoided.componentA, avoided.componentB))

    assert avoided_pair not in resulting_pairs
