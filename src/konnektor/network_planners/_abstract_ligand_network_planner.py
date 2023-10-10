import abc
import logging
import itertools
from typing import Iterable, Tuple, List

from openfe.setup import LigandNetwork

from gufe import SmallMoleculeComponent, AtomMapper

log = logging.getLogger(__name__)
#log.setLevel(logging.WARNING)

class easyLigandNetworkPlanner(abc.ABC):
    def __init__(self, mapper, scorer, network_generator):
        self.mapper = mapper
        self.scorer =  scorer
        self.network_generator = network_generator

    def __call__(self, *args, **kwargs):
        return self.generate_ligand_network(*args, **kwargs)

    def _input_generate_all_possible_mappings(self, ligands)->Tuple:
        #Todo: this should actually become the function of max network.

        # First create a network with all the proposed mappings (scored)
        mapping_generator = itertools.chain.from_iterable(
            self.mapper.suggest_mappings(molA, molB)
            for molA, molB in itertools.combinations(ligands, 2))
        mappings = [mapping.with_annotations({'score': self.scorer(mapping)})
                    for mapping in mapping_generator]

        #mapping_scores = [mapping.annotations['score'] for mapping in mappings]
        return ligands, mappings,


    def generate_ligand_network(self, ligands)->LigandNetwork:
        """Plan a Network which connects all ligands with minimal cost

        Parameters
        ----------
        ligands : Iterable[SmallMoleculeComponent]
        the ligands to include in the Network
        """
        nodes, edges, weights = self._input_generate_all_possible_mappings(ligands=ligands)
        selected_edges = self.network_generator(nodes=nodes, edges=edges, weights=weights)

        return LigandNetwork(edges=selected_edges, nodes=ligands)


    def concatenate_networks(self, ligandNetworks:List[LigandNetwork], nEdges:int = 3) -> LigandNetwork:

        selected_edges = []
        selected_nodes = []
        connecting_edges = []

        log.info("Number of edges in individual networks:\n"+str(sum([len(s.edges) for s in ligandNetworks]))+"/"+str([len(s.edges) for s in ligandNetworks]))
        for ligandNetworkA, ligandNetworkB in itertools.combinations(ligandNetworks, 2):
            compoundsA = list(ligandNetworkA.nodes)
            compoundsB = list(ligandNetworkB.nodes)

            #bipartite Graph
            mappings = []
            for cA in compoundsA:
                for cB in compoundsB:
                    m = [mapping.with_annotations({'score': self.scorer(mapping)}) for mapping in self.mapper.suggest_mappings(cA,cB)]
                    mappings.extend(m)

            #prio Queue
            mappings = list(sorted(mappings, key=lambda x: x.annotations["score"]))
            connecting_nodes = []

            #nEdges diversity increase:
            for mapping in mappings:
                nodes = [mapping.componentA.name, mapping.componentB.name]
                if(all([not n in connecting_nodes for n in nodes])):
                    connecting_edges.append(mapping)
                    connecting_nodes.extend(nodes)
                    if(len(connecting_edges)>=nEdges):
                        break

            log.info("Adding ConnectingEdges: "+str(len(connecting_edges)))

            selected_edges.extend(ligandNetworkA.edges)
            selected_edges.extend(ligandNetworkB.edges)
            selected_nodes.extend(compoundsA)
            selected_nodes.extend(compoundsB)

        selected_edges = list(set(selected_edges))
        selected_edges.extend(connecting_edges)
        log.info("Total Concatenating Edges: "+str(len(connecting_edges)))
        log.info("Total Concatenated Edges: "+str(len(selected_edges)))

        concat_LigandNetwork = LigandNetwork(edges=selected_edges, nodes=set(selected_nodes))
        return concat_LigandNetwork
