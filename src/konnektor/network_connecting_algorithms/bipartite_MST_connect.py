import logging
import itertools
from typing import List
from konnektor.utils import LigandNetwork


log = logging.getLogger(__name__)

class mst_concatenate():
    def __init__(self, mapper, scorer):
        self.mapper = mapper
        self.scorer = scorer
    def concatenate_networks(self, ligand_networks: List[LigandNetwork], n_edges: int = 3) -> LigandNetwork:

        selected_edges = []
        selected_nodes = []
        connecting_edges = []

        log.info("Number of edges in individual networks:\n " +str(sum([len(s.edges) for s in ligand_networks]) ) +"/ " +str
            ([len(s.edges) for s in ligand_networks]))
        for ligandNetworkA, ligandNetworkB in itertools.combinations(ligand_networks, 2):
            compoundsA = list(ligandNetworkA.nodes)
            compoundsB = list(ligandNetworkB.nodes)

            # bipartite Graph
            mappings = []
            for cA in compoundsA:
                for cB in compoundsB:
                    m = [mapping.with_annotations({'score': self.scorer(mapping)}) for mapping in self.mapper.suggest_mappings(cA ,cB)]
                    mappings.extend(m)

            # prio Queue
            mappings = list(sorted(mappings, key=lambda x: x.annotations["score"]))
            connecting_nodes = []

            # n_edges diversity increase:
            for mapping in mappings:
                nodes = [mapping.componentA.name, mapping.componentB.name]
                if(all([not n in connecting_nodes for n in nodes])):
                    connecting_edges.append(mapping)
                    connecting_nodes.extend(nodes)
                    if(len(connecting_edges ) >= n_edges):
                        break

            log.info("Adding ConnectingEdges:  " +str(len(connecting_edges)))

            selected_edges.extend(ligandNetworkA.edges)
            selected_edges.extend(ligandNetworkB.edges)
            selected_nodes.extend(compoundsA)
            selected_nodes.extend(compoundsB)

        selected_edges = list(set(selected_edges))
        selected_edges.extend(connecting_edges)
        log.info("Total Concatenating Edges:  " +str(len(connecting_edges)))
        log.info("Total Concatenated Edges:  " +str(len(selected_edges)))

        concat_LigandNetwork = LigandNetwork(edges=selected_edges, nodes=set(selected_nodes))
        return concat_LigandNetwork
