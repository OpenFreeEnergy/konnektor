from konnektor.network_tools.merge import merge_two_networks, merge


from konnektor.utils.toy_data import build_random_mst_network, build_two_random_mst_network
from konnektor.network_analysis import get_is_connected


def test_merge_two_mst_networks():
    networkA, networkB = build_two_random_mst_network(n_compounds=20, sub_networks=2, rand_seed=42)

    new_network = merge_two_networks(networkA, networkB)

    assert len(new_network.nodes) == 20
    assert len(new_network.edges) == 19
    assert get_is_connected(new_network)


def test_merge_mst_networks():
    networks = build_two_random_mst_network(n_compounds=20, sub_networks=2, rand_seed=42)

    new_network = merge(networks)

    assert len(networks) == 2
    assert len(new_network.nodes) == 20
    assert len(new_network.edges) == 19
    assert get_is_connected(new_network)

    networks = build_two_random_mst_network(n_compounds=20, sub_networks=4, overlap=1, rand_seed=42)
    new_network = merge(networks)

    assert len(networks) == 4
    assert len(new_network.nodes) == 20
    assert len(new_network.edges) == 19
    assert get_is_connected(new_network)