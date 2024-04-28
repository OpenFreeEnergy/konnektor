import itertools
import logging
from typing import Iterable

from konnektor.network_planners.generators.netx_netgen import MstNetworkGenerator

from ._abstract_ligand_network_concatenator import LigandNetworkConcatenator
from ...utils import LigandNetwork

log = logging.getLogger(__name__)


# Todo: check this algorithm again

class MstConcatenate(LigandNetworkConcatenator):
    def __init__(self, mapper: AtomMapper, scorer: AtomMappingScorer, n_connecting_edges: int = 3, nprocesses: int = 1):
        """
        This concatenator is connnecting two Networks with a kruskal like approach up to the number of connecting edges.

        Parameters
        ----------
        mapper: AtomMapper
            the atom mapper is required, to define the connection between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        n_connecting_edges: int, optional
            number of connecting edges. (default: 3)
        nprocesses: int
            number of processes that can be used for the network generation. (default: 1)
        """
        super().__init__(mapper=mapper, scorer=scorer, network_generator=MstNetworkGenerator(), nprocesses=nprocesses)
        self.n_connecting_edges = n_connecting_edges

    def concatenate_networks(self, ligand_networks: Iterable[LigandNetwork]) -> LigandNetwork:
        """

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            an iterable of ligand networks, that shall be connected.
        n_connecting_edges: int
            number of edges, to connect the networks

        Returns
        -------
        LigandNetwork
            returns a concatenated LigandNetwork object, containing all networks.

        """

        log.info("Number of edges in individual networks:\n " + str(sum([len(s.edges) for s in ligand_networks])) +
                 "/ " + str([len(s.edges) for s in ligand_networks]))
        selected_edges = []
        selected_nodes = []
        for ligandNetworkA, ligandNetworkB in itertools.combinations(ligand_networks, 2):
            # Generate Bipartite Graph
            # Todo: use max network here!
            ligands = []
            bipartite_graph_mappings = []
            for cA in list(ligandNetworkA.nodes):
                for cB in list(ligandNetworkB.nodes):
                    m = [mapping.with_annotations({'score': self.scorer(mapping), "type": "connect"})
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
            log.info("Adding ConnectingEdges:  " + str(len(selected_mappings)))

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

            # Constructed final Edges:
            # Add all old network edges:
            for network in ligand_networks:
                selected_edges.extend(network.edges)
                selected_nodes.extend(network.nodes)

            # Add network connecting edges
            selected_edges.extend(selected_mappings)

        concat_LigandNetwork = LigandNetwork(edges=selected_edges, nodes=set(selected_nodes))

        log.info("Total Concatenated Edges:  " + str(len(selected_edges)))

        return concat_LigandNetwork
