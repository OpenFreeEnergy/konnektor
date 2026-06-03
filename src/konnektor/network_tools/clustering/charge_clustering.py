# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from collections import defaultdict

from gufe import SmallMoleculeComponent

from ._abstract_clusterer import _AbstractClusterer


class ChargeClusterer(_AbstractClusterer):
    """Clusters molecules based on their formal charges."""

    def __init__(self):
        pass

    def cluster_compounds(
        self,
        components: list[SmallMoleculeComponent],
    ) -> dict[int, list[SmallMoleculeComponent]]:
        """Cluster compounds according to formal charge.

        Parameters
        ----------
        components : list[SmallMoleculeComponent]

        Returns
        -------
        dict[int, list[gufe.SmallMoleculeComponent]]
            Dict of formal charge states mapped to lists of the corresponding molecules.
        """
        clusters = defaultdict(list)

        for m in components:
            rdk_m = m.to_rdkit()
            fc = sum(a.GetFormalCharge() for a in rdk_m.GetAtoms())

            clusters[fc].append(m)

        return dict(clusters)
