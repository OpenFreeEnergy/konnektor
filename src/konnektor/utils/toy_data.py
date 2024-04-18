from rdkit import Chem
from rdkit.Chem import AllChem
from gufe import SmallMoleculeComponent

from gufe import LigandAtomMapping
import numpy as np


def build_random_dataset(n=20):
    class genMapper:
        def suggest_mappings(self, molA, molB):
            yield LigandAtomMapping(molA, molB, {})

    class genScorer:
        vals = np.random.random(int(n ** 2))
        i = 0

        def __call__(self, x):
            v = self.vals[self.i]
            self.i = (self.i + 1) % n
            return v

    mols = [Chem.AddHs(Chem.MolFromSmiles("C" * i)) for i in range(1, int(n))]
    [Chem.rdDistGeom.EmbedMolecule(m) for m in mols]
    compounds = [SmallMoleculeComponent(name=str(i), rdkit=m) for i, m in enumerate(mols)]

    return compounds, genMapper, genScorer
