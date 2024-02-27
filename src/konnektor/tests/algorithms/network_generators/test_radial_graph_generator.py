# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
import numpy as np
import networkx as nx

from konnektor.network_planners.generators.netx_netgen import RadialNetworkGenerator
from konnektor.tests.data.conf import nine_mols_edges


def test_radial_network_generation_find_center(nine_mols_edges):

        edges = [(e[0], e[1]) for e in nine_mols_edges]
        weights = [e[2] for e in nine_mols_edges]

        gen = RadialNetworkGenerator()
        c_node, avg_weight = gen._central_lig_selection(edges, weights)

        assert c_node == "lig_14" # Check central node
        assert np.round(avg_weight,2) == 2.86

def test_radial_network_generation_without_center(nine_mols_edges):

        edges = [(e[0], e[1]) for e in nine_mols_edges]
        weights = [e[2] for e in nine_mols_edges]
        nodes = set([n for e in edges for n in e])

        gen = RadialNetworkGenerator()
        g = gen.generate_network(edges, weights)

        assert len(nodes)-1 == len(g.edges)
        assert all(["lig_14" in e for e in g.edges]) #check central node
        assert all([e[0]!=e[1] for e in g.edges]) # No self connectivity
        assert isinstance(g, nx.Graph)


def test_radial_network_generation_with_center(nine_mols_edges):

        edges = [(e[0], e[1]) for e in nine_mols_edges]
        weights = [e[2] for e in nine_mols_edges]
        nodes = set([n for e in edges for n in e])

        gen = RadialNetworkGenerator()
        g = gen.generate_network(edges, weights, central_node="lig_11")

        assert len(nodes)-1 == len(g.edges)
        assert all(["lig_11" in e for e in g.edges]) #check central node
        assert all([e[0]!=e[1] for e in g.edges]) # No self connectivity
        assert isinstance(g, nx.Graph)

