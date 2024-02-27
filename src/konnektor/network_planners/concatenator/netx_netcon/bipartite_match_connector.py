import logging
import networkx as nx


log = logging.getLogger(__name__)

class MatchingConcatenator():

    def concatenate_networks(self, nodesA, nodesB,
                             edges, weights)->nx.Graph:
        """
        TODO Separate networking from Ligand stuff
        Parameters
        ----------
        ligand_networks
        n_connecting_edges

        Returns
        -------

        """

        wedges_map = {(e[0], e[1]): w for e,w in zip(edges, weights)}
        wedges = [(e[0], e[1], w) for e,w in zip(edges, weights)]

        #Build Bi-partite graph
        B = nx.Graph()
        B.add_nodes_from(nodesA, bipartite=0)
        B.add_nodes_from(nodesB, bipartite=1)
        B.add_weighted_edges_from(wedges)

        #Matching
        m_edges = [(e[0], e[1], wedges_map[e]) for e in nx.maximal_matching(B)]
        mg = nx.Graph()
        mg.add_nodes_from(nodesA, bipartite=0)
        mg.add_nodes_from(nodesB, bipartite=1)
        mg.add_weighted_edges_from(m_edges)

        return mg
