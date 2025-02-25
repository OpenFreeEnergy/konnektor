# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor
import itertools

import networkx as nx

from konnektor.network_planners._networkx_implementations import (
    bipartite_match_algorithm,
)
from konnektor.tests.data.conf import nine_mols_edges


def test_mst_network_generation(nine_mols_edges):
    nodes = list(set([n for e in nine_mols_edges for n in e]))
    nodesA = nodes[:4]
    nodesB = nodes[4:]
    edges = list(itertools.product(nodesA, nodesB))
    weights = [1 for i in range(len(edges))]

    gen = bipartite_match_algorithm.MatchingConcatenatAlgorithm()
    g = gen.concatenate_networks(nodesA, nodesB, edges, weights)

    print(g.edges)
    assert isinstance(g, nx.Graph)
    assert len(nodesA) == len(g.edges)
    assert all([e[0] != e[1] for e in g.edges])  # No self connectivity
