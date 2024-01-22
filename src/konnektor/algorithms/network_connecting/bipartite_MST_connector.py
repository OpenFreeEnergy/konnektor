import logging
import itertools
from typing import List
from konnektor.utils import LigandNetwork


log = logging.getLogger(__name__)

class MstConcatenate():
    def __init__(self, mapper, scorer,  n_connecting_edges: int = 3):
        self.mapper = mapper
        self.scorer = scorer
        self.n_connecting_edges = n_connecting_edges
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
        selected_edges = []
        selected_nodes = []
        connecting_edges = []

        log.info("Number of edges in individual networks:\n " + str(sum([len(s.edges) for s in ligand_networks])) +
                 "/ " + str([len(s.edges) for s in ligand_networks]))
        for ligandNetworkA, ligandNetworkB in itertools.combinations(ligand_networks, 2):

            # Generate Bipartite Graph
            bipartite_graph_edges = []
            for cA in list(ligandNetworkA.nodes):
                for cB in list(ligandNetworkB.nodes):
                    m = [mapping.with_annotations({'score': self.scorer(mapping), "type":"connect"})
                         for mapping in self.mapper.suggest_mappings(cA, cB)]
                    bipartite_graph_edges.extend(m)

            # prio Queue -
            connecting_nodes = []
            priority_queue = list(sorted(bipartite_graph_edges, key=lambda x: x.annotations["score"]))
            for mapping in priority_queue:
                nodes = [mapping.componentA.name, mapping.componentB.name]

                # increase connecting node diversity:
                if all([n not in connecting_nodes for n in nodes]):
                    setattr(mapping, "type", "connect") #give edge flavor
                    connecting_edges.append(mapping)
                    connecting_nodes.extend(nodes)

                    if len(connecting_edges) >= self.n_connecting_edges:
                        break

            log.info("Adding ConnectingEdges:  " + str(len(connecting_edges)))

        # Constructed final Edges:
        # Add all old network edges:
        for network in ligand_networks:
            selected_edges.extend(network.edges)
            selected_nodes.extend(network.nodes)

        # Add network connecting edges
        selected_edges.extend(connecting_edges)

        concat_LigandNetwork = LigandNetwork(edges=selected_edges, nodes=set(selected_nodes))

        log.info("Total Concatenating Edges:  " + str(len(connecting_edges)))
        log.info("Total Concatenated Edges:  " + str(len(selected_edges)))

        return concat_LigandNetwork
