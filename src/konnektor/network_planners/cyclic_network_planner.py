from typing import Union, List, Iterable

from gufe import SmallMoleculeComponent
from konnektor.utils import LigandNetwork

from ..network_generator_algorithms import CyclicNetworkGenerator
from ._abstract_ligand_network_planner import LigandNetworkPlanner
from .maximal_network_planner import MaximalNetworkPlanner


class CyclicLigandNetworkPlanner(LigandNetworkPlanner):

    def __init__(self, mapper, scorer, node_present_in_cycles=2, cycle_sizes:Union[int, List[int]]=3 ):
        super().__init__(mapper=mapper, scorer=scorer,
                       network_generator=CyclicNetworkGenerator(node_cycle_connectivity=node_present_in_cycles,
                                                              sub_cycle_size_range=cycle_sizes),
                         _initial_edge_lister=MaximalNetworkPlanner(
                             mapper=mapper, scorer=scorer))


    def generate_ligand_network(self, ligands: Iterable[SmallMoleculeComponent]) ->LigandNetwork:
            # Build Full Graph
            initial_networks = self._initial_edge_lister.generate_ligand_network(
                nodes=ligands)
            mappings = initial_networks.edges

            # Translate Mappings to graphable:
            print("prepare network")
            edge_map = {(ligands.index(m.componentA), ligands.index(m.componentB)): m for m in mappings}
            edges = list(sorted(edge_map.keys()))
            weights = [edge_map[k].annotations['score'] for k in edges]

            print("calculate Network")
            cg = self.network_generator.generate_network_double_greedy(edges=edges, weights=weights)

            selected_mappings = [edge_map[k] for k in cg.edges]
            print("Done")

            return LigandNetwork(edges=selected_mappings, nodes=ligands)
