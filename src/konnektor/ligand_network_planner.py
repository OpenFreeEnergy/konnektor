from typing import Iterable

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