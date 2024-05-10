from collections import defaultdict

import gufe

from ._abstract_clusterer import _AbstractClusterer


class ChargeClusterer(_AbstractClusterer):
    """Clusters molecules based on their formal charges"""

    def __init__(self):
        pass

    def cluster_compounds(self,
                          components: list[gufe.SmallMoleculeComponent]) -> \
    dict[int, list[gufe.SmallMoleculeComponent]]:
        """Cluster compounds according to formal charge

        Returns a dict which has keys of all present formal charge states,
        mapping to lists of the corresponding molecules
        """
        clusters = defaultdict(list)

        for m in components:
            rdk_m = m.to_rdkit()
            fc = sum(a.GetFormalCharge() for a in rdk_m.GetAtoms())

            clusters[fc].append(m)

        return dict(clusters)
