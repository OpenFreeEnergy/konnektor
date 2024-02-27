import itertools
from typing import Iterable

from gufe import SmallMoleculeComponent
from konnektor.utils import LigandNetwork

from .netx_netgen import RadialNetworkGenerator
from ._abstract_ligand_network_planner import LigandNetworkPlanner
from .maximal_network_planner import MaximalNetworkPlanner

class StarLigandNetworkPlanner(LigandNetworkPlanner):

    def __init__(self, mapper, scorer):
        super().__init__(mapper=mapper, scorer=scorer, network_generator=RadialNetworkGenerator(),
                         _initial_edge_lister=MaximalNetworkPlanner(
                             mapper=mapper, scorer=scorer))


    def generate_ligand_network(self, ligands: Iterable[SmallMoleculeComponent],
                                central_ligand:SmallMoleculeComponent=None) -> LigandNetwork:
            # Build Full Graph
            ligands = list(ligands)

            if(central_ligand is None):
                #Full Graph Construction
                initial_network = self._initial_edge_lister.generate_ligand_network(
                    nodes=ligands)
                mappings = initial_network.edges

                #Translate Mappings to graphable:
                edge_map = {(ligands.index(m.componentA), ligands.index(m.componentB)): m for m in mappings}
                edges = list(sorted(edge_map.keys()))
                weights = [edge_map[k].annotations['score'] for k in edges]

                rg = self.network_generator.generate_network(edges=edges, weights=weights)
                selected_mappings = [edge_map[k] for k in rg.edges]

            else:   #Given central ligands: less effort. - Trivial Case

                if self.scorer is None:
                    scorer = lambda x: -1
                else:
                    scorer = self.scorer

                mapping_generators = [self.mapper.suggest_mappings(
                    central_ligand, molA) for molA in ligands]
                selected_mappings = [mapping.with_annotations({'score': scorer(mapping)})
                            for mapping_generator in mapping_generators for
                                     mapping in mapping_generator]

            return LigandNetwork(edges=selected_mappings, nodes=ligands)

RadialLigandNetworkPlanner = StarLigandNetworkPlanner