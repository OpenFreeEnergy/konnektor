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
        self.nBits = 1

    @property
    def fpSize(self):
        return self.nBits

    # Scikit-Learn expects to be able to set fpSize directly on object via .set_params(), so this updates nBits used by the abstract class
    @fpSize.setter
    def fpSize(self, fpSize):
        self.nBits = fpSize

    def _mol2fp(self, mol):
        pass

    def _transform_mol(self, mol):
        return np.array([Chem.GetFormalCharge(mol)])
