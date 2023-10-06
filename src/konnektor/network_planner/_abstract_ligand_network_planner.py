import abc
import itertools
from typing import Iterable, Tuple

from openfe.setup import LigandNetwork

from gufe import SmallMoleculeComponent

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


'''
class ligandNetworkPlanner():
    def __int__(self, mappers=[], scorers=[], network_planner=None):
        self.mappers = mappers
        self.mapping_scorer = scorers
        self.network_planner = network_planner

    def __call__(self, *args, **kwargs):
        self.generate_ligand_network()

    def generate_ligand_network(self, ligands:Iterable):

        mappings = [[[mapper(molA, molB) for molB in ligands[i:]] for i, molA in enumerate(ligands)] for mapper in self.mappers]
        scores = [[[(mapping_scorer(mapping) for mapping_scorer in self.mapping_scorer) for mapping in molA_mappings] for molA_mappings in mapper_mappings] for mapper_mappings in mappings]

        #flatten_stuff
        mappings_scores =[]
        for mapping_mappper, mapping_scores in zip(mappings, scores):
            for mol_mappings, mol_score in zip(mapping_mappper, mapping_scores):
                mappings_scores.extend(list(zimp(mol_mappings, mol_score)))

        self.network_planner()

'''