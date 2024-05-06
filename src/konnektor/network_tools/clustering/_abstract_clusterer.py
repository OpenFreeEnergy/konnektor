import abc
from collections import defaultdict
import gufe


class _AbstractClusterer:
    """General interface for clustering"""

    @abc.abstractmethod
    def cluster_compounds(self, components: list[gufe.SmallMoleculeComponent]) -> dict[int, list[gufe.SmallMoleculeComponent]]:
        """Cluster compounds according to a given algorithm """
