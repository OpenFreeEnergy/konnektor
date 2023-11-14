from typing import Iterable

from gufe import SmallMoleculeComponent
from konnektor.utils import LigandNetwork

from src.konnektor.algorithms.network_generator_algorithms import RadialNetworkGenerator
from ._abstract_ligand_network_planner import easyLigandNetworkPlanner


class RadialLigandNetworkPlanner(easyLigandNetworkPlanner):

    def __init__(self, mapper, scorer):
        super().__init__(mapper=mapper, scorer=scorer, network_generator=RadialNetworkGenerator())



    def generate_ligand_network(self, ligands: Iterable[SmallMoleculeComponent],
                                central_ligand:SmallMoleculeComponent=None) -> LigandNetwork:
            # Build Full Graph
            ligands = list(ligands)

            if(central_ligand is None):
                #Full Graph Construction
                ligands, mappings = self._input_generate_all_possible_mappings(ligands=ligands)

                #Translate Mappings to graphable:
                edge_map = {(ligands.index(m.componentA), ligands.index(m.componentB)): m for m in mappings}
                edges = list(sorted(edge_map.keys()))
                weights = [edge_map[k].annotations['score'] for k in edges]

                rg = self.network_generator.generate_network(edges=edges, weights=weights)
                selected_mappings = [edge_map[k] for k in rg.edges]

            else:   #Given central ligands: less effort. - Trivial Case
                mapping_generator = [self.mapper.suggest_mappings(central_ligand, molA) for molA in ligands]
                selected_mappings = [mapping.with_annotations({'score': self.scorer(mapping)})
                            for mapping in mapping_generator]

            return LigandNetwork(edges=selected_mappings, nodes=ligands)
