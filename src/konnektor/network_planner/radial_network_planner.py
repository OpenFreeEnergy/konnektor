# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/kartograf

import math
import numpy as np
import networkx as nx
from typing import Iterable, Callable
import itertools

from gufe import AtomMapper, AtomMapping
from gufe import SmallMoleculeComponent

from openfe.setup import LigandNetwork    # only temproary
from ._abstract_network_planner import _AbstractNetworkPlanner, Network


class radialNetworkPlanner(_AbstractNetworkPlanner):

    def __init__(self,  metric_aggregation_method:Callable=None):
        self.metric_aggregation_method = metric_aggregation_method


    def _central_lig_selection(self, edges, weights)->int:
        nodes = set([n for e in edges for n in e])
        edge_weights = list(zip(edges, weights))

        node_scores = {n: [e_s[1] for e_s in edge_weights if (n in e_s[0])] for n in nodes}
        aggregated_scores = list(map(lambda x: (x[0], np.sum(x[1])), node_scores.items()))
        sorted_node_scores = list(sorted(aggregated_scores, key=lambda x: x[1]))

        opt_node = sorted_node_scores[0]
        return opt_node

    def generate_network(self,edges, weights, central_node=None) -> nx.Graph:
        """Generate a radial network with all ligands connected to a central node

        Also known as hub and spoke or star-map, this plans a Network where
        all ligands are connected via a central ligand.

        Parameters
        ----------
        ligands : iterable of SmallMoleculeComponents
          the ligands to arrange around the central ligand
        central_ligand : SmallMoleculeComponent
          the ligand to use as the hub/central ligand
        mappers : iterable of AtomMappers
          mappers to use, at least 1 required
        scorer : scoring function, optional
          a callable which returns a float for any AtomMapping.  Used to
          assign scores to potential mappings, higher scores indicate worse
          mappings.

        Raises
        ------
        ValueError
          if no mapping between the central ligand and any other ligand can be
          found

        Returns
        -------
        network : Network
          will have an edge between each ligand and the central ligand, with the
          mapping being the best possible mapping found using the supplied atom
          mappers.
          If no scorer is supplied, the first mapping provided by the iterable
          of mappers will be used.
        """

        if(central_node is None):
            central_node = self._central_lig_selection(edges=edges, weights=weights)

        radial_edges = []
        for edge, weight in zip(edges):
            if(central_node in edge):
                radial_edges.append([edge[0], edge[1], weight])

        rg = nx.Graph()
        rg.add_weighted_edges_from(ebunch_to_add=edges)

        return rg
