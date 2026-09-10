# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
from gufe import LigandNetwork

from konnektor.network_planners.concatenators import MstConcatenator
from konnektor.utils.toy_data import (
    EmptyMapper,
    RandomScorer,
    build_n_random_mst_network,
    build_random_dataset,
)


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


def test_score_inter_network_edges_respects_exclusions():
    """The exclude_edges should never get scored."""
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
    mappings = concatenator._score_inter_network_edges(networkA, networkB, exclude=set())
    excluded = frozenset((mappings[0].componentA, mappings[0].componentB))

    # Re-run with that mapping excluded.
    mappings = concatenator._score_inter_network_edges(networkA, networkB, exclude={excluded})
    resulting_pairs = {frozenset((mapping.componentA, mapping.componentB)) for mapping in mappings}

    assert excluded not in resulting_pairs


def test_tied_scores_pick_highest_key():
    """On a score tie, the mapping with the max .key is chosen (deterministic)."""
    n = 20
    networkA, networkB = build_n_random_mst_network(
        n_compounds=n, sub_networks=2, overlap=0, rand_seed=42
    )

    def constant_scorer(mapping):
        return 0.5

    concatenator = MstConcatenator(EmptyMapper(), constant_scorer)

    # every candidate between the two subnetworks scores 0.5
    candidates = concatenator._score_inter_network_edges(networkA, networkB, exclude=set())
    expected = max(candidates, key=lambda m: m.key)

    bridges = concatenator._select_mst_bridges([networkA, networkB], exclude=set())
    assert len(bridges) == 1
    assert bridges[0] == expected


def test_spanning_tree_pairs_discards_partial_forest():
    """A forest that does not span all subnetworks is discarded."""
    n = 4
    components, _, _ = build_random_dataset(n_compounds=n, rand_seed=1)
    subnetworks = [LigandNetwork(nodes=[c], edges=[]) for c in components]
    concatenator = MstConcatenator(EmptyMapper(), RandomScorer(n=n))

    # Only subnetworks 0, 1, and 2 are connected by candidate mappings.
    mapping_01 = concatenator._score_inter_network_edges(subnetworks[0], subnetworks[1], exclude=[])[0]
    mapping_12 = concatenator._score_inter_network_edges(subnetworks[1], subnetworks[2], exclude=[])[0]
    best_mapping_by_pair = {(0, 1): mapping_01, (1, 2): mapping_12}

    # Two edges are sufficient to span three subnetworks, but not four.
    partial = concatenator._select_spanning_tree_pairs(best_mapping_by_pair, n_networks=4)
    assert partial == []
    # Sanity check: the same candidate graph is a complete tree over three.
    complete = concatenator._select_spanning_tree_pairs(best_mapping_by_pair, n_networks=3)
    assert len(complete) == 2
