# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from typing import Callable, Iterable

import networkx as nx
import numpy as np
from gufe import SmallMoleculeComponent

from ._abstract_network_algorithm import _AbstractNetworkAlgorithm


class RadialNetworkAlgorithm(_AbstractNetworkAlgorithm):
    def __init__(self, metric_aggregation_method: Callable = None, n_centers: int = 1):
        self.metric_aggregation_method = metric_aggregation_method
        self.n_centers = n_centers

    def _central_lig_selection(
        self, edges: list[tuple[int, int]], weights: list[float]
    ) -> Iterable[int]:
        nodes = set([n for e in edges for n in e])
        # The initial "weights" are Scores, which need to be translated to weights.
        weights = list(map(lambda x: 1 - x, weights))
        edge_weights = list(zip(edges, weights))

        node_scores = {
            n: [e_s[1] for e_s in edge_weights if (n in e_s[0])] for n in nodes
        }
        filtered_node_scores = dict(
            filter(lambda x: len(x[1]) != len(nodes), node_scores.items())
        )

        if len(filtered_node_scores) == 0:  # Todo: allow relaxed criterion here.
            raise ValueError("Could not find a single node connecting all edges!")

        aggregated_scores = list(
            map(lambda x: (x[0], np.sum(x[1])), filtered_node_scores.items())
        )
        sorted_node_scores = list(sorted(aggregated_scores, key=lambda x: x[1]))
        opt_nodes = sorted_node_scores[: self.n_centers]
        return opt_nodes

    def generate_network(
        self,
        edges: list[tuple[int, int]],
        weights: list[float],
        central_node: int = None,
    ) -> nx.Graph:
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

        if central_node is None:
            central_nodes = self._central_lig_selection(
                edges=edges,
                weights=weights,
            )
        elif isinstance(central_node, (SmallMoleculeComponent, str)):
            central_nodes = [(central_node, 1)]
        else:
            raise ValueError("invalide central node type: " + str(type(central_node)))

        wedges = []
        for edge, weight in zip(edges, weights):
            if any(central_node in edge for central_node, avg_score in central_nodes):
                wedges.append([edge[0], edge[1], weight])

        # Todo: Warning if something was not connected to the central ligand?
        nodes = set([n for e in edges for n in e])

        self.radial_graph = nx.Graph()
        [self.radial_graph.add_node(n) for n in nodes]
        self.radial_graph.add_weighted_edges_from(
            ebunch_to_add=[(e[0], e[1], e[2]) for e in wedges]
        )

        return self.radial_graph
