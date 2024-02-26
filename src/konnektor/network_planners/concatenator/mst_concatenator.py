import logging
import itertools
from typing import List
from konnektor.utils import LigandNetwork

from ...algorithms.network_generator import MstNetworkGenerator

log = logging.getLogger(__name__)

class MstConcatenate():
    def __init__(self, mapper, scorer,  n_connecting_edges: int = 3):
        self.mapper = mapper
        self.scorer = scorer
        self.n_connecting_edges = n_connecting_edges
        self.network_generator =MstNetworkGenerator()
    def concatenate_networks(self, ligand_networks: List[LigandNetwork]) -> LigandNetwork:
        """
        TODO Separate networking from Ligand stuff
        Parameters
        ----------
        ligand_networks
        n_connecting_edges

        Returns
        -------

        """

        log.info("Number of edges in individual networks:\n " + str(sum([len(s.edges) for s in ligand_networks])) +
                 "/ " + str([len(s.edges) for s in ligand_networks]))
        selected_edges = []
        selected_nodes = []
        for ligandNetworkA, ligandNetworkB in itertools.combinations(ligand_networks, 2):
            # Generate Bipartite Graph
            ligands = []
            bipartite_graph_mappings = []
            for cA in list(ligandNetworkA.nodes):
                for cB in list(ligandNetworkB.nodes):
                    m = [mapping.with_annotations({'score': self.scorer(mapping), "type":"connect"})
                         for mapping in self.mapper.suggest_mappings(cA, cB)]
                    bipartite_graph_mappings.extend(m)
                    ligands.extend([cA, cB])

            # Find MST subset for Pipartite
            edge_map = {
                (ligands.index(m.componentA), ligands.index(m.componentB)): m
                for m in bipartite_graph_mappings}
            edges = list(edge_map.keys())
            weights = [edge_map[k].annotations['score'] for k in edges]

            mg = self.network_generator.generate_network(edges, weights,
                                                         n_edges=self.n_connecting_edges)

            selected_mappings = [edge_map[k] if (k in edge_map) else edge_map[
                    tuple(list(k)[::-1])] for k in mg.edges]

            """ # prio queue kruska approach
            connecting_nodes = []
            priority_queue = list(sorted(bipartite_graph_mappings, key=lambda x: x.annotations["score"]))
            for mapping in priority_queue:
                nodes = [mapping.componentA.name, mapping.componentB.name]
    
                # increase connecting node diversity:
                if all([n not in connecting_nodes for n in nodes]):
                    setattr(mapping, "type", "connect") #give edge flavor
                    connecting_edges.append(mapping)
                    connecting_nodes.extend(nodes)
    
                    if len(connecting_edges) >= self.n_connecting_edges:
                        break
            """
            log.info("Adding ConnectingEdges:  " + str(len(selected_mappings)))

            # Constructed final Edges:
            # Add all old network edges:
            for network in ligand_networks:
                selected_edges.extend(network.edges)
                selected_nodes.extend(network.nodes)

            # Add network connecting edges
            selected_edges.extend(selected_mappings)

        concat_LigandNetwork = LigandNetwork(edges=selected_edges, nodes=set(selected_nodes))

        log.info("Total Concatenating Edges:  " + str(len(selected_mappings)))
        log.info("Total Concatenated Edges:  " + str(len(selected_edges)))

        return concat_LigandNetwork
