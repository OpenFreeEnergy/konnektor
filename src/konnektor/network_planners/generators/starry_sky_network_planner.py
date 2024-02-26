from typing import Union, List, Iterable

from gufe import SmallMoleculeComponent
from konnektor.utils import LigandNetwork

from ...algorithms.network_generator import StarrySkyNetworkGenerator
from ._abstract_ligand_network_planner import LigandNetworkPlanner
from .maximal_network_planner import MaximalNetworkPlanner

class StarrySkyLigandNetworkPlanner(LigandNetworkPlanner):

    def __init__(self, mapper, scorer, target_node_connectivity: int = 3):
        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=StarrySkyNetworkGenerator(target_node_connectivity=target_node_connectivity),
                         _initial_edge_lister=MaximalNetworkPlanner(mapper=mapper, scorer=scorer))

    def generate_ligand_network(self, ligands: Iterable[SmallMoleculeComponent]) -> LigandNetwork:
        # Build Full Graph
        initial_networks = self._initial_edge_lister.generate_ligand_network(
            nodes=ligands)
        mappings = initial_networks.edges

        # Translate Mappings to graphable:
        edge_map = {(ligands.index(m.componentA), ligands.index(m.componentB)): m for m in mappings}
        edges = list(sorted(edge_map.keys()))
        weights = [edge_map[k].annotations['score'] for k in edges]

        sg = self.network_generator.generate_network(edges=edges, weights=weights)

        selected_mappings = [edge_map[k] if(k in edge_map) else edge_map[tuple(list(k)[::-1])] for k in sg.edges]
        return LigandNetwork(edges=selected_mappings, nodes=ligands)
