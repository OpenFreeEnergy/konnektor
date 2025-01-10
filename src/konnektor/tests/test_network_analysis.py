# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
import numpy as np

from gufe import LigandNetwork
from konnektor.network_analysis import (
    get_network_score,
    get_is_connected,
    get_component_connectivities,
    get_component_number_cycles,
    get_number_of_network_cycles,
    get_component_scores,
    get_transformation_failure_robustness,
)

from konnektor.utils.toy_data import (
    build_random_mst_network,
    build_random_fully_connected_network,
)


#   Graph Topology related functionalities - Connectivity
def test_get_is_connected():
    # Check for connected.
    g = build_random_mst_network()
    assert get_is_connected(g)

    # Check for disconnected.
    g_disscon = LigandNetwork(nodes=g.nodes, edges=list(g.edges)[1:])
    assert get_is_connected(g_disscon) == False


def test_get_node_connectives():
    n_compounds = 30
    expected_set = set(map(str, range(n_compounds)))
    expected_arr = np.zeros(shape=n_compounds)
    expected_arr += n_compounds - 1

    g = build_random_fully_connected_network(n_compounds)

    cons = get_component_connectivities(g, normalize=False)

    assert len(cons) == n_compounds
    assert len(set(cons.keys()).intersection(expected_set)) == 0
    assert all(n == (len(g.nodes) - 1) for i, n in cons.items())
    np.testing.assert_array_almost_equal(expected_arr, [n for i, n in cons.items()])

    cons = get_component_connectivities(g, normalize=True)
    expected_arr = np.zeros(shape=n_compounds)
    expected_arr += (n_compounds - 1) / (n_compounds * (n_compounds - 1) / 2)

    assert len(cons) == n_compounds
    assert len(set(cons.keys()).intersection(expected_set)) == 0
    np.testing.assert_array_almost_equal(expected_arr, [n for i, n in cons.items()])


@pytest.mark.flaky(reruns=3)  # pytest-rerunfailures;
@pytest.mark.parametrize(
    "grapher,expected_robustness,failure_rate,exact",
    [
        (
            build_random_fully_connected_network,
            1,
            0.05,
            True,
        ),  # never gets disconnected
        (
            build_random_fully_connected_network,
            0,
            1,
            True,
        ),  # all edges removed, should always fail
        (build_random_mst_network, 0.0, 0.05, True),  # Mst always fails
    ],
)
def test_get_edge_failure_robustness(grapher, expected_robustness, failure_rate, exact):
    n_compounds = 5
    g = grapher(n_compounds, rand_seed=42)

    robustness = get_transformation_failure_robustness(
        g, failure_rate=failure_rate, nrepeats=100, seed=42
    )

    if exact:
        np.testing.assert_allclose(actual=robustness, desired=expected_robustness)
    else:
        np.testing.assert_allclose(actual=robustness, desired=expected_robustness, rtol=0.25)


#   Graph Topology related functionalities - Cycles
@pytest.mark.parametrize("n_cycle_size", [3, 4])
def test_get_node_number_cycles_fully_connected_graph(n_cycle_size):
    n_compounds = 30
    expected_set = set(map(str, range(n_compounds)))

    g = build_random_fully_connected_network(n_compounds)

    # each node has $\sum^(node_cycle_size)_(i) n_nodes-i$ (for cycle size 3 ==> 29 * 28) different possible node combinatoins for cycles.
    expected_cycles_per_cs = []
    for cycle_size in range(3, n_cycle_size + 1):
        expected_ncyles = np.prod([n_compounds - i for i in range(1, cycle_size)]) / 2
        expected_cycles_per_cs.append(expected_ncyles)
    expected_number_of_cycles_per_node = sum(expected_cycles_per_cs)

    result_dict = get_component_number_cycles(g, higher_bound=n_cycle_size)

    assert len(result_dict) == n_compounds
    assert len(set(result_dict.keys()).intersection(expected_set)) == 0

    expected_arr = np.zeros(shape=n_compounds)
    expected_arr += expected_number_of_cycles_per_node
    np.testing.assert_array_almost_equal(x=expected_arr, y=[n for i, n in result_dict.items()])


def test_get_node_number_cycles_mst_graph():
    n_compounds = 30
    n_cycle_size = 3
    expected_set = set(map(str, range(n_compounds)))
    g = build_random_mst_network(n_compounds)

    # should not contain cycles as all
    result_dict = get_component_number_cycles(g, higher_bound=n_cycle_size)

    assert len(result_dict) == n_compounds
    assert len(set(result_dict.keys()).intersection(expected_set)) == 0

    expected_arr = np.zeros(shape=n_compounds)  # expect zero cycles in a tree
    np.testing.assert_array_almost_equal(expected_arr, [n for i, n in result_dict.items()])


@pytest.mark.parametrize("n_cycle_size", [3, 4])
def test_get_number_of_graph_cycles_fully_connected_graph(n_cycle_size):
    n_compounds = 30
    g = build_random_fully_connected_network(n_compounds)

    # Compute the number of cycles as before, except, now make each cycle unique (if cycle size == 3 make sure to count only 1)
    expected_cycles_per_cs = []
    for cycle_size in range(3, n_cycle_size + 1):
        expected_ncyles = n_compounds * np.prod([n_compounds - i for i in range(1, cycle_size)]) / 2
        expected_cycles_per_cs.append(expected_ncyles / cycle_size)
    expected_number_of_graph_cycles = sum(expected_cycles_per_cs)

    number_of_graph_cycles = get_number_of_network_cycles(g, higher_bound=n_cycle_size)

    assert expected_number_of_graph_cycles == number_of_graph_cycles


@pytest.mark.parametrize("n_cycle_size", [3, 4])
def test_get_number_of_graph_cycles_mst_graph(n_cycle_size):
    n_compounds = 30
    g = build_random_mst_network(n_compounds)

    # no cycle should occur
    expected_number_of_graph_cycles = 0

    number_of_graph_cycles = get_number_of_network_cycles(g, higher_bound=n_cycle_size)

    assert expected_number_of_graph_cycles == number_of_graph_cycles


# Graph Edge Score related Functionality - convoluted scores


def test_get_mst_graph_score():
    # Check for graph scores.
    seed = 42
    n_compounds = 30
    g = build_random_mst_network(n_compounds=n_compounds, rand_seed=seed)
    np.testing.assert_allclose(get_network_score(g), 27.52, atol=1e-3)


def test_get_fully_connected_graph_score():
    # Check for graph scores.
    g = build_random_fully_connected_network()
    np.testing.assert_allclose(get_network_score(g), 191.629, atol=1e-3)


def test_get_norm_node_scores_fully_connected_graph():
    n_compounds = 30
    expected_set = set(map(str, range(n_compounds)))
    g = build_random_fully_connected_network(n_compounds, uni_score=True)

    n_scores = get_component_scores(g, normalize=True)

    expected_arr = np.zeros(shape=n_compounds)
    expected_arr += (n_compounds - 1) / ((n_compounds - 1) * n_compounds) * 2

    assert len(n_scores) == n_compounds
    assert len(set(n_scores.keys()).intersection(expected_set)) == 0
    np.testing.assert_array_almost_equal(x=expected_arr, y=[n for i, n in n_scores.items()])


def test_get_node_scores_fully_connected_graph():
    n_compounds = 30
    expected_set = set(map(str, range(n_compounds)))
    g = build_random_fully_connected_network(n_compounds, uni_score=True)

    n_scores = get_component_scores(g, normalize=False)

    expected_arr = np.zeros(shape=n_compounds)
    expected_arr += n_compounds - 1

    assert len(n_scores) == n_compounds
    assert len(set(n_scores.keys()).intersection(expected_set)) == 0
    np.testing.assert_array_almost_equal(x=expected_arr, y=[n for i, n in n_scores.items()])


def test_get_node_scores_fully_connected_graph():
    n_compounds = 30
    expected_set = set(map(str, range(n_compounds)))
    g = build_random_fully_connected_network(n_compounds, uni_score=True)

    n_scores = get_component_scores(g, normalize=False)

    expected_arr = np.zeros(shape=n_compounds)
    expected_arr += n_compounds - 1

    assert len(n_scores) == n_compounds
    assert len(set(n_scores.keys()).intersection(expected_set)) == 0
    np.testing.assert_array_almost_equal(x=expected_arr, y=[n for i, n in n_scores.items()])
