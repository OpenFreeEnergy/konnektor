# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import itertools
import networkx as nx

from typing import Iterable, Callable

from gufe import AtomMapper, AtomMapping
from gufe import SmallMoleculeComponent

from openfe.setup.ligand_network import LigandNetwork    # only temproary
from ._abstract_network_generator import _AbstractNetworkGenerator, Network

class MstNetworkGenerator(_AbstractNetworkGenerator):


    def generate_network(self, edges, weights) -> nx.Graph:
        wedges = []
        for edge, weight in zip(edges, weights):
            wedges.append([edge[0], edge[1], weight])

        self.g = g = nx.Graph()
        g.add_weighted_edges_from(ebunch_to_add=wedges)

        # Next analyze that network to create minimal spanning network. Because
        # we carry the original (directed) AtomMapping, we don't lose
        # direction information when converting to an undirected graph.
        min_edges = nx.minimum_spanning_edges(g, weight='score')

        mse = [edge_data for edge_data in min_edges]

        mg = nx.Graph()
        mg.add_weighted_edges_from(ebunch_to_add=mse)

        return mg
