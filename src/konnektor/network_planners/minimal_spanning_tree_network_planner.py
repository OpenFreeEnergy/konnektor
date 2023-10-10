from openfe import LigandNetwork

from network_generator_algorithms import MstNetworkGenerator
from network_planners._abstract_ligand_network_planner import easyLigandNetworkPlanner


class MinimalSpanningTreeLigandNetworkPlanner(easyLigandNetworkPlanner):

    def __init__(self, mapper, scorer):
        """Plan a Network which connects all ligands with minimal cost

        Parameters
        ----------
        mappers : Iterable[AtomMapper]
        the AtomMappers to use to propose mappings.  At least 1 required,
        but many can be given, in which case all will be tried to find the
        lowest score edges
        scorer : Scoring function
        any callable which takes a AtomMapping and returns a float
        """
        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=MstNetworkGenerator())

    def generate_ligand_network(self, ligands) ->LigandNetwork:


        ligands, mappings = self._input_generate_all_possible_mappings(ligands=ligands)

        # Translate Mappings to graphable:
        edge_map = {(ligands.index(m.componentA), ligands.index(m.componentB)): m for m in mappings}
        edges = list(edge_map.keys())
        weights = [edge_map[k].annotations['score'] for k in edges]

        mg = self.network_generator.generate_network(edges, weights)

        if(len(mg.nodes) < len(ligands)):
            nodes_index = {l:ligands.index(l) for l in ligands}
            missing_nodes = [l for l in ligands if(nodes_index[l] in mg.nodes)]
            raise RuntimeError("Unable to create edges for some nodes: "
                               + str(list(missing_nodes)))

        selected_mappings = [edge_map[k] if(k in edge_map) else edge_map[tuple(list(k)[::-1])] for k in mg.edges]

        return LigandNetwork(edges=selected_mappings, nodes=ligands)
