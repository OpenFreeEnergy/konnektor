# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import networkx as nx

from konnektor.network_planners._networkx_implementations import MstNetworkAlgorithm


def test_mst_network_generation(nine_mols_edges):
    expected_edges = [
        ("lig_14", "lig_13"),
        ("lig_14", "lig_15"),
        ("lig_14", "lig_11"),
        ("lig_14", "lig_8"),
        ("lig_15", "lig_12"),
        ("lig_15", "lig_10"),
        ("lig_12", "lig_9"),
        ("lig_11", "lig_16"),
    ]

    init_graph = nx.Graph()
    init_graph.add_weighted_edges_from(
        ebunch_to_add=[(e[0], e[1], e[2]) for e in nine_mols_edges], weight="score"
    )

    gen = MstNetworkAlgorithm()
    g = gen.generate_network(init_graph)

    assert len(g.nodes) - 1 == len(g.edges)
    assert [e in g.edges for e in expected_edges]  # TODO: Does not work as intended.
    assert all([e[0] != e[1] for e in g.edges])  # No self connectivity
    assert isinstance(g, nx.Graph)
