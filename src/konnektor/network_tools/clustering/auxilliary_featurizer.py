# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from typing import Union

import numpy as np
from rdkit import Chem
from scikit_mol.fingerprints import FpsTransformer


class ChargeTransformer(FpsTransformer):
    def __init__(self, parallel: Union[bool, int] = False):
        """Calculates the RDKit FormalCharge and provides it as single field vector.

        Parameters
        ----------
        parallel: Union[bool, int], optional
        Whether to parallelize the calculations, default: False

        """
        super().__init__(parallel=parallel)
        self.fpSize = 1

    def _mol2fp(self, mol):
        pass

    def _transform_mol(self, mol):
        return np.array([Chem.GetFormalCharge(mol)])
