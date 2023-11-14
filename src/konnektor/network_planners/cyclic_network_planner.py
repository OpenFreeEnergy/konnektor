from typing import Union, List, Iterable

from gufe import SmallMoleculeComponent
from konnektor.utils import LigandNetwork

from src.konnektor.algorithms.network_generator_algorithms import CyclicNetworkGenerator
from ._abstract_ligand_network_planner import easyLigandNetworkPlanner


class CyclicLigandNetworkPlanner(easyLigandNetworkPlanner):

    def __init__(self, mapper, scorer, node_present_in_cycles=2, cycle_sizes:Union[int, List[int]]=3 ):
        super().__init__(mapper=mapper, scorer=scorer,
                       network_generator=CyclicNetworkGenerator(node_cycle_connectivity=node_present_in_cycles,
                                                              sub_cycle_size_range=cycle_sizes))


    def generate_ligand_network(self, ligands: Iterable[SmallMoleculeComponent]) ->LigandNetwork:
            # Build Full Graph
            ligands, mappings = self._input_generate_all_possible_mappings(ligands=ligands)

            # Translate Mappings to graphable:
            edge_map = {(ligands.index(m.componentA), ligands.index(m.componentB)): m for m in mappings}
            edges = list(sorted(edge_map.keys()))
            weights = [edge_map[k].annotations['score'] for k in edges]

            cg = self.network_generator.generate_network(edges=edges, weights=weights)

            selected_mappings = [edge_map[k] for k in cg.edges]
            return LigandNetwork(edges=selected_mappings, nodes=ligands)
