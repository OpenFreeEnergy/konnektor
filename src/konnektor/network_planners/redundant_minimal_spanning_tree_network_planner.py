from  konnektor.utils import LigandNetwork

from ..algorithms.network_generator import MstNetworkGenerator
from ._abstract_ligand_network_planner import LigandNetworkPlanner
from .maximal_network_planner import MaximalNetworkPlanner


class RedundantMinimalSpanningTreeLigandNetworkPlanner(LigandNetworkPlanner):

    def __init__(self, mapper, scorer, n_redundancy:int=3):
        """Plan a Network which connects all ligands n times with minimal cost.
        This planner uses n_redundancy times the MST algorithm on the full
        graph, disallowing to use an already selected edge twice for the MST calculation.

        Parameters
        ----------
        mappers : Iterable[AtomMapper]
            the AtomMappers to use to propose mappings.  At least 1 required,
            but many can be given, in which case all will be tried to find the
            lowest score edges
        scorer : Scoring function
            any callable which takes a AtomMapping and returns a float
        n_redundancy: int
            use MST n times to get a redundant set of edges.

        """
        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=MstNetworkGenerator(),
                         _initial_edge_lister=MaximalNetworkPlanner(
                             mapper=mapper, scorer=scorer))

        self.n_redundancy = n_redundancy

    def generate_ligand_network(self, ligands) ->LigandNetwork:


        initial_network = self._initial_edge_lister.generate_ligand_network(
            nodes=ligands)
        mappings = initial_network.edges

        # Translate Mappings to graphable:
        edge_map = {(ligands.index(m.componentA), ligands.index(m.componentB)): m for m in mappings}
        edges = list(edge_map.keys())
        weights = [edge_map[k].annotations['score'] for k in edges]

        edge_weight = dict(zip(edges, weights))

        selected_edges=[]
        for _ in range(self.n_redundancy):
            edges = list(edge_weight.keys())

            # filter for already selected edges
            filter_sEdges = lambda x: (x not in selected_edges and
                                       x[::-1] not in selected_edges)
            edges = list(filter(filter_sEdges, edges))
            #print("Edges", len(edges), edges)

            weights = [edge_weight[e] for e in edges]
            edge_weight = dict(list(zip(edges, weights)))
            #print("dEdges", len(edge_weight))
            ns = set([n for e in edges for n in e])
            #print("nodes", len(ns), ns)

            mg = self.network_generator.generate_network(edges, weights)

            if not mg.connected:
                nodes_index = {l:ligands.index(l) for l in ligands}
                missing_nodes = [l for l in ligands if(nodes_index[l] in mg.nodes)]
                raise RuntimeError("Unable to create edges for some nodes: "
                                   + str(list(missing_nodes)))
            #print("sel", len(mg.edges), mg.edges)

            selected_edges.extend(list(mg.edges))
            #print("collected", len(selected_edges), selected_edges)

        selected_mappings = [edge_map[k] if(k in edge_map) else edge_map[
            tuple(list(k)[::-1])] for k in selected_edges]

        return LigandNetwork(edges=selected_mappings, nodes=ligands)
