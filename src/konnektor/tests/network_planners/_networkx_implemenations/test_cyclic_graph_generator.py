# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import networkx as nx

from konnektor.network_planners._networkx_implementations import \
    CyclicNetworkAlgorithm
from konnektor.tests.data.conf import nine_mols_edges


# more test here also for the params
def test_cyclic_network_generation(nine_mols_edges):
    expected_edges = [('lig_15', 'lig_12'), ('lig_15', 'lig_16'),
                      ('lig_15', 'lig_13'),
                      ('lig_15', 'lig_14'), ('lig_15', 'lig_10'),
                      ('lig_15', 'lig_11'),
                      ('lig_15', 'lig_9'), ('lig_15', 'lig_8'),
                      ('lig_8', 'lig_13'),
                      ('lig_8', 'lig_14'), ('lig_9', 'lig_14'),
                      ('lig_9', 'lig_13'),
                      ('lig_16', 'lig_13'), ('lig_16', 'lig_14'),
                      ('lig_12', 'lig_13'),
                      ('lig_12', 'lig_14'), ('lig_13', 'lig_14'),
                      ('lig_13', 'lig_11'),
                      ('lig_13', 'lig_10'), ('lig_10', 'lig_14'),
                      ('lig_14', 'lig_11')]

    edges = [(e[0], e[1]) for e in nine_mols_edges]
    weights = [e[2] for e in nine_mols_edges]
    nodes = set([n for e in edges for n in e])

    gen = CyclicNetworkAlgorithm(sub_cycle_size_range=[3],
                                 node_cycle_connectivity=2)
    g = gen.generate_network(edges, weights)

    print(g.edges)
    assert (len(nodes) - 1) * 2 < len(g.edges)  # min size
    assert all([len([e for e in g.edges if (n in e)]) >= 2 for n in
                nodes])  # minimal required edges
    assert [e in g.edges for e in expected_edges]
    assert all([e[0] != e[1] for e in g.edges])  # No self connectivity
    assert isinstance(g, nx.Graph)

def test_cyclic_dg_network_generation(nine_mols_edges):
    expected_edges = [('lig_15', 'lig_12'), ('lig_15', 'lig_16'),
                      ('lig_15', 'lig_13'),
                      ('lig_15', 'lig_14'), ('lig_15', 'lig_10'),
                      ('lig_15', 'lig_11'),
                      ('lig_15', 'lig_9'), ('lig_15', 'lig_8'),
                      ('lig_8', 'lig_13'),
                      ('lig_8', 'lig_14'), ('lig_9', 'lig_14'),
                      ('lig_9', 'lig_13'),
                      ('lig_16', 'lig_13'), ('lig_16', 'lig_14'),
                      ('lig_12', 'lig_13'),
                      ('lig_12', 'lig_14'), ('lig_13', 'lig_14'),
                      ('lig_13', 'lig_11'),
                      ('lig_13', 'lig_10'), ('lig_10', 'lig_14'),
                      ('lig_14', 'lig_11')]

    edges = [(e[0], e[1]) for e in nine_mols_edges]
    weights = [e[2] for e in nine_mols_edges]
    nodes = set([n for e in edges for n in e])

    gen = CyclicNetworkAlgorithm(sub_cycle_size_range=[3],
                                 node_cycle_connectivity=2)
    g = gen.generate_network_double_greedy(edges, weights)

    print(g.edges)
    assert (len(nodes) - 1) * 2 < len(g.edges)  # min size
    assert all([len([e for e in g.edges if (n in e)]) >= 2 for n in
                nodes])  # minimal required edges
    assert [e in g.edges for e in expected_edges]
    assert all([e[0] != e[1] for e in g.edges])  # No self connectivity
    assert isinstance(g, nx.Graph)
