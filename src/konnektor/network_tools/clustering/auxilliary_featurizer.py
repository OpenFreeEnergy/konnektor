# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor


import numpy as np
from rdkit import Chem
from scikit_mol.fingerprints.baseclasses import BaseFpsTransformer


class ChargeTransformer(BaseFpsTransformer):
    def __init__(self, n_jobs: int | None = None):
        """Calculates the RDKit FormalCharge and provides it as single field vector.

        Parameters
        ----------
        n_jobs : int, optional default=None
            The maximum number of concurrently running jobs.
            None is a marker for 'unset' that will be interpreted as n_jobs=1 unless the call is performed under a parallel_config()
        Whether to parallelize the calculations, default: False

        """
        super().__init__(n_jobs=n_jobs)
        self.fpSize = 1

    def _transform_mol(self, mol):
        return np.array([Chem.GetFormalCharge(mol)])
