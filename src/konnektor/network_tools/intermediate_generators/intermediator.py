
from rdkit import Chem

# Intermediate generator
from generator import generation
from generator import scoring
#notes
# Clean up imports
# rename package.
# how reliant on openeye?
# Deappreceation Morgan Generator

from gufe import SmallMoleculeComponent


class Intermediator():
    def __init__(self, scoring_method=None, fp_type:str="ECFP4", num_tries:int=5,
                 num_random_smiles: int= 10, collect_bidirectional:bool=True,
                 exponent_path:int=4, n_rounds:int=1, num_random_samples:int=500,
                 num_mutation_ls=(1,2),
                 exponent_local_chemical_space:int= 2,contribution_lomap:float= 0.2, contribution_similarity:float= 0.8):
        self.scoring_method = scoring_method
        self.num_tries = num_tries
        self.num_random_smiles = num_random_smiles
        self.collect_bidirectional = collect_bidirectional
        self.exponent_path = exponent_path
        self.fp_type=fp_type
        self.n_rounds = n_rounds

        self.num_random_samples =num_random_samples
        self.num_mutation_ls =num_mutation_ls

        self.exponent_local_chemical_space = exponent_local_chemical_space
        self.contribution_lomap = contribution_lomap
        self.contribution_similarity= contribution_similarity

    def generate_intermediate(self, molA: SmallMoleculeComponent, molB: SmallMoleculeComponent):

        rdmolA = molA.to_rdkit()
        rdmolB = molB.to_rdkit()
        liga_smiles = Chem.MolToSmiles(rdmolA)
        ligb_smiles = Chem.MolToSmiles(rdmolB)

        print("> Phase I")
        if self.scoring_method == "3D":
            generated_paths = generation.generate_multiple_paths_rocs(
                liga_smiles, ligb_smiles, self.num_tries, self.num_random_smiles,
                collect_bidirectional=self.collect_bidirectional,
                exponent_path=self.exponent_path, n_rounds=self.n_rounds,
                fp_type=self.fp_type)
        else:
            generated_paths = generation.generate_multiple_paths(liga_smiles,
                                                                 ligb_smiles,
                                                                 self.num_tries,
                                                                 self.num_random_smiles,
                                                                 collect_bidirectional=self.collect_bidirectional,
                                                                 exponent_path=self.exponent_path,
                                                                 n_rounds=self.n_rounds,
                                                                 fp_type=self.fp_type)
            print(generated_paths)
        print("> Phase II")

        generated_mols = generation.generate_chemical_space(liga_smiles,
                                                            ligb_smiles,
                                                            generated_paths,
                                                            self.num_random_samples,
                                                            self.num_mutation_ls,
                                                            fp_type=self.fp_type)
        print("> Phase III")

        if self.scoring_method == "3D":
            sorted_smiles_dict = scoring.score_molecules_lomap_rocs(liga_smiles,
                                                                    ligb_smiles,
                                                                    generated_mols,
                                                                    self.exponent_local_chemical_space,
                                                                    self.contribution_lomap,
                                                                    self.contribution_similarity)
        else:
            sorted_smiles_dict = scoring.score_molecules_lomap_tanimoto(
                liga_smiles, ligb_smiles, generated_mols,
                self.exponent_local_chemical_space, self.contribution_lomap,
                self.contribution_similarity)

        print(sorted_smiles_dict)
        print("> Phase IV")

        selected_intermediate = next(iter(sorted_smiles_dict))
        rd_mol_intermediate = Chem.MolFromSmiles(selected_intermediate)
        return SmallMoleculeComponent.from_rdkit(rdkit=rd_mol_intermediate,
                                                 name=molA.name+"_"+molB.name+"_intermediate")