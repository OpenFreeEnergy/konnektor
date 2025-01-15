import pytest

from konnektor.network_analysis import get_is_connected
from konnektor.network_tools.network_handling.delete import (
    delete_component,
    delete_transformation,
)
from konnektor.utils.toy_data import (
    build_random_fully_connected_network,
    build_random_mst_network,
)


def test_delete_fc_component():
    network = build_random_fully_connected_network(n_compounds=10)
    del_node = list(network.nodes)[0]

    new_network = delete_component(network=network, component=del_node)

    assert len(new_network.nodes) == len(network.nodes) - 1
    assert len(new_network.edges) < len(network.edges)
    assert del_node not in new_network.nodes
    assert get_is_connected(new_network)


def test_delete_connected_mst_component():
    network = build_random_mst_network(n_compounds=10)
    # get an connecting edges
    edges = [(e.componentA, e.componentB) for e in list(network.edges)]
    n_count = {n: sum(1 if n in e else 0 for e in edges) for n in list(network.nodes)}
    del_node = None

    for n, c in n_count.items():
        if c > 1:
            del_node = n
            break

    with pytest.raises(RuntimeError, match="Resulting network is not connected anymore!"):
        new_network = delete_component(network=network, component=del_node)


def test_delete_mst_component():
    network = build_random_mst_network(n_compounds=10)
    del_node = list(network.nodes)[0]

    new_network = delete_component(network=network, component=del_node, must_stay_connected=False)

    assert len(new_network.nodes) == len(network.nodes) - 1
    assert len(new_network.edges) < len(network.edges)
    assert del_node not in new_network.nodes


def test_delete_fc_transformation():
    network = build_random_fully_connected_network(n_compounds=10)
    del_edge = list(network.edges)[0]

    new_network = delete_transformation(network=network, transformation=del_edge)

    assert len(new_network.nodes) == len(network.nodes)
    assert len(new_network.edges) == len(network.edges) - 1
    assert del_edge not in new_network.edges
    assert get_is_connected(new_network)


def test_delete_connected_mst_transformation():
    network = build_random_mst_network(n_compounds=10)
    del_edge = list(network.edges)[0]

    with pytest.raises(RuntimeError, match="Resulting network is not connected anymore!"):
        new_network = delete_transformation(network=network, transformation=del_edge)


def test_delete_mst_transformation():
    network = build_random_mst_network(n_compounds=10)
    del_edge = list(network.edges)[0]

    new_network = delete_transformation(
        network=network, transformation=del_edge, must_stay_connected=False
    )

    assert len(new_network.nodes) == len(network.nodes)
    assert len(new_network.edges) == len(network.edges) - 1
    assert del_edge not in new_network.edges
    assert not get_is_connected(new_network)
