# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import abc

import gufe


class _AbstractClusterer:
    """General interface for clustering"""

    @abc.abstractmethod
    def cluster_compounds(
        self, components: list[gufe.SmallMoleculeComponent]
    ) -> dict[int, list[gufe.SmallMoleculeComponent]]:
        """Cluster compounds according to a given algorithm"""
