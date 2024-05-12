# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import networkx as nx

from konnektor.network_planners._networkx_implementations import \
    NNodeEdgesNetworkAlgorithm
from konnektor.tests.data.conf import nine_mols_edges


def test_n_edge_network_generation(nine_mols_edges):
    expected_edges = [('lig_14', 'lig_13'), ('lig_14', 'lig_15'),
                      ('lig_14', 'lig_11'), ('lig_14', 'lig_8'),
                      ('lig_15', 'lig_12'), ('lig_15', 'lig_10'),
                      ('lig_12', 'lig_9'), ('lig_11', 'lig_16')]

    edges = [(e[0], e[1]) for e in nine_mols_edges]
    weights = [e[2] for e in nine_mols_edges]
    nodes = set([n for e in edges for n in e])

    n_edges = 2
    gen = NNodeEdgesNetworkAlgorithm(target_node_connectivity=n_edges)
    g = gen.generate_network(edges, weights)

    assert isinstance(g, nx.Graph)
    assert len(nodes) * n_edges >= len(g.edges)
    assert [e in g.edges for e in expected_edges]
    assert all([e[0] != e[1] for e in g.edges])  # No self connectivity
