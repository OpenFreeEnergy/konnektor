# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
import networkx as nx
import numpy as np

from konnektor.network_planners._networkx_implementations import RadialNetworkAlgorithm


def test_radial_network_generation_find_center(nine_mols_edges):
    edges = [(e[0], e[1]) for e in nine_mols_edges]
    weights = [e[2] for e in nine_mols_edges]

    gen = RadialNetworkAlgorithm()
    c_node, avg_weight = gen._central_lig_selection(edges, weights)[0]

    assert c_node == "lig_10"  # Check central node
    np.testing.assert_allclose(avg_weight, 2.055, rtol=0.01)


@pytest.mark.parametrize("n_centers", [2, 3, 4])
def test_radial_network_generation_find_centers(nine_mols_edges, n_centers):
    edges = [(e[0], e[1]) for e in nine_mols_edges]
    weights = [e[2] for e in nine_mols_edges]

    gen = RadialNetworkAlgorithm(n_centers=n_centers)
    centers = gen._central_lig_selection(edges, weights)

    print(centers)
    expected_centers = ["lig_10", "lig_8", "lig_9", "lig_16"]
    expected_weights = [
        2.0551095953189917,
        3.6524873109359146,
        4.270400420741822,
        4.543886935944357,
    ]
    for i, (cID, avg_weight) in enumerate(centers):
        print(cID, avg_weight)
        assert cID == expected_centers[i]  # Check central node
        np.testing.assert_allclose(avg_weight, expected_weights[i], rtol=0.01)


def test_radial_network_generation_without_center(nine_mols_edges):
    edges = [(e[0], e[1]) for e in nine_mols_edges]
    weights = [e[2] for e in nine_mols_edges]
    nodes = set([n for e in edges for n in e])

    gen = RadialNetworkAlgorithm()
    g = gen.generate_network(edges, weights)

    assert len(nodes) - 1 == len(g.edges)
    assert all(["lig_10" in e for e in g.edges])  # check central node
    assert all([e[0] != e[1] for e in g.edges])  # No self connectivity
    assert isinstance(g, nx.Graph)


def test_radial_network_generation_with_center(nine_mols_edges):
    edges = [(e[0], e[1]) for e in nine_mols_edges]
    weights = [e[2] for e in nine_mols_edges]
    nodes = set([n for e in edges for n in e])

    gen = RadialNetworkAlgorithm()
    g = gen.generate_network(edges, weights, central_node="lig_11")

    assert len(nodes) - 1 == len(g.edges)
    assert all(["lig_11" in e for e in g.edges])  # check central node
    assert all([e[0] != e[1] for e in g.edges])  # No self connectivity
    assert isinstance(g, nx.Graph)
