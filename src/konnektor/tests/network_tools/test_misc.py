from konnektor.network_analysis import get_is_connected
from konnektor.network_tools.misc import delete_component, delete_transformation
from konnektor.utils.toy_data import build_random_mst_network, \
    build_random_fully_connected_network


def test_delete_fc_component():
    network = build_random_fully_connected_network(n_compounds=10)
    del_node = list(network.nodes)[0]

    new_network = delete_component(network=network, component=del_node)

    assert len(new_network.nodes) == len(network.nodes) - 1
    assert len(new_network.edges) < len(network.edges)
    assert del_node not in new_network.nodes
    assert get_is_connected(new_network)


def test_delete_mst_component():
    network = build_random_mst_network(n_compounds=30)
    del_node = list(network.nodes)[0]

    new_network = delete_component(network=network, component=del_node)

    assert len(new_network.nodes) == len(network.nodes) - 1
    assert len(new_network.edges) < len(network.edges)
    assert del_node not in new_network.nodes
    assert not get_is_connected(new_network)


def test_delete_fc_transformation():
    network = build_random_fully_connected_network(n_compounds=10)
    del_edge = list(network.edges)[0]

    new_network = delete_transformation(network=network,
                                        transformation=del_edge)

    assert len(new_network.nodes) == len(network.nodes)
    assert len(new_network.edges) == len(network.edges) - 1
    assert del_edge not in new_network.edges
    assert get_is_connected(new_network)


def test_delete_mst_transformation():
    network = build_random_mst_network(n_compounds=10)
    del_edge = list(network.edges)[0]

    new_network = delete_transformation(network=network,
                                        transformation=del_edge)

    assert len(new_network.nodes) == len(network.nodes)
    assert len(new_network.edges) == len(network.edges) - 1
    assert del_edge not in new_network.edges
    assert not get_is_connected(new_network)
