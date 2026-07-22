# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import abc

from gufe import SmallMoleculeComponent


class _AbstractClusterer:
    """General interface for clustering"""

    @abc.abstractmethod
    def cluster_compounds(
        self, components: list[SmallMoleculeComponent]
    ) -> dict[int, list[SmallMoleculeComponent]]:
        """Cluster compounds according to a given algorithm"""
