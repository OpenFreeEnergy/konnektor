# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from gufe import LigandNetwork

from konnektor.network_tools.network_handling.decompose import decompose_network
from konnektor.utils.toy_data import build_n_random_mst_network


def test_decompose_connected_network_returns_single():
    network = build_n_random_mst_network(n_compounds=20, rand_seed=42)
    assert network.is_connected()

    sub_networks = decompose_network(network)

    assert len(sub_networks) == 1
    assert sub_networks[0].nodes == network.nodes
    assert sub_networks[0].edges == network.edges


def test_decompose_disconnected_network():
    networkA, networkB = build_n_random_mst_network(
        n_compounds=20,
        sub_networks=2,
        overlap=0,
        rand_seed=42,
    )
    networks = [networkA, networkB]
    nodes = set().union(*(n.nodes for n in networks))
    edges = set().union(*(n.edges for n in networks))
    disconnected = LigandNetwork(nodes=nodes, edges=edges)
    assert not disconnected.is_connected()

    sub_networks = decompose_network(disconnected)

    assert len(sub_networks) == 2
    assert all(sn.is_connected() for sn in sub_networks)
    # nothing is lost or duplicated
    assert set().union(*(sn.nodes for sn in sub_networks)) == disconnected.nodes
    assert set().union(*(sn.edges for sn in sub_networks)) == disconnected.edges
    assert sum(len(sn.nodes) for sn in sub_networks) == len(disconnected.nodes)
