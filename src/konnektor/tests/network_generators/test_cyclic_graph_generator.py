# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
import numpy as np
import networkx as nx

from konnektor.network_generator_algorithms import CyclicNetworkGenerator
from konnektor.tests.data.conf import nine_mols_edges


def test_mst_network_generation_find_center(nine_mols_edges):
        expected_edges = [('lig_15', 'lig_12'), ('lig_15', 'lig_16'), ('lig_15', 'lig_13'),
                          ('lig_15', 'lig_14'), ('lig_15', 'lig_10'), ('lig_15', 'lig_11'),
                          ('lig_15', 'lig_9'), ('lig_15', 'lig_8'), ('lig_8', 'lig_13'),
                          ('lig_8', 'lig_14'), ('lig_9', 'lig_14'), ('lig_9', 'lig_13'),
                          ('lig_16', 'lig_13'), ('lig_16', 'lig_14'), ('lig_12', 'lig_13'),
                          ('lig_12', 'lig_14'), ('lig_13', 'lig_14'), ('lig_13', 'lig_11'),
                          ('lig_13', 'lig_10'), ('lig_10', 'lig_14'), ('lig_14', 'lig_11')]

        edges = [(e[0], e[1]) for e in nine_mols_edges]
        weights = [e[2] for e in nine_mols_edges]
        nodes = set([n for e in edges for n in e])

        gen = CyclicNetworkGenerator(sub_cycle_size_range=[3], node_cycle_connectivity=2)
        g = gen.generate_network(edges, weights)

        print(g.edges)
        assert (len(nodes)-1)*2 < len(g.edges) #min size
        assert all([len([e for e in g.edges if(n in e)])>=2 for n in nodes]) #minimal required edges
        assert [e in g.edges for e in expected_edges]
        assert all([e[0] != e[1] for e in g.edges]) # No self connectivity
        assert isinstance(g, nx.Graph)
