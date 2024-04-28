from typing import Iterable
from rdkit import Chem
from rdkit.Chem import AllChem
from gufe import SmallMoleculeComponent

from gufe import LigandAtomMapping, AtomMapper, AtomMappingScorer, AtomMapping
import numpy as np


class genMapper(AtomMapper):
    def __init__(self):
        pass

    def suggest_mappings(self, molA, molB)->AtomMapping:
        yield LigandAtomMapping(molA, molB, {})

class genScorer(AtomMappingScorer):
    def __init__(self, n_scores:int):
        vals = np.random.random(int(n_scores ** 2))
        i = 0

    def __call__(self, x)->float:
        v = self.vals[self.i]
        self.i = (self.i + 1) % self.n_scores
        return v


def build_random_dataset(n=20)-> (Iterable[SmallMoleculeComponent], AtomMapper, AtomMappingScorer):

    gen_mapper = genMapper()
    gen_scorer = genScorer(n_scores=n)

    #generate random Molecules
    smiles = ["".join(np.random.choice(["C", "O", "N", "S"], replace=True, size= i%30)) for i in range(1, int(n))]
    mols = [Chem.AddHs(Chem.MolFromSmiles(smiles)) for i in smiles]
    [Chem.rdDistGeom.EmbedMolecule(m) for m in mols]
    compounds = [SmallMoleculeComponent(name=str(i), rdkit=m) for i, m in enumerate(mols)]

    return compounds, gen_mapper, gen_scorer
