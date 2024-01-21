from typing import Iterable, Callable, Union, Optional

from gufe import SmallMoleculeComponent, AtomMapper
from konnektor.utils import LigandNetwork    # only temproary

from ._abstract_ligand_network_planner import LigandNetworkPlanner

from sklearn.cluster import KMeans
from scikit_mol.fingerprints import RDKitFingerprintTransformer, MorganFingerprintTransformer
from scikit_mol.descriptors import MolecularDescriptorTransformer

from konnektor.utils import CompoundDiversityClustering
from konnektor.network_connecting_algorithms.bipartite_MST_connect import mst_concatenate

class DiversityNetworkPlanner(LigandNetworkPlanner):
    def __init__(self, mapper, scorer,
                 node_featurizer=RDKitFingerprintTransformer(),
                 clustering= KMeans(n_clusters=3)):
        super().__init__(mapper=mapper, scorer=scorer,
                       network_generator=None)
        self.feat =  node_featurizer
        self.cluster =clustering

    def generate_ligand_network(self,
                         nodes: Iterable[SmallMoleculeComponent],
                         progress: Union[bool, Callable[[Iterable], Iterable]] = True,
                         # allow_disconnected=True
                         ):
        from konnektor.network_planners import MinimalSpanningTreeLigandNetworkPlanner, CyclicLigandNetworkPlanner

        #step 1: Seperate nodes by diversity
        cc = CompoundDiversityClustering(featurize=self.feat,
                                         cluster=self.cluster)

        clusters = cc.cluster_compounds(nodes)


        # Sub Network, based on clusters
        planner = CyclicLigandNetworkPlanner(mapper=self.mapper, scorer=self.scorer)
        alt_planner = MinimalSpanningTreeLigandNetworkPlanner(mapper=self.mapper, scorer=self.scorer)

        sub_networks = []
        for cID, mols in clusters.items():
            if(len(mols)>1):
                if(len(mols)>2):
                    sub_network = planner(mols)
                    sub_networks.append(sub_network)
                else:
                    sub_network = alt_planner(mols)
                    sub_networks.append(sub_network)
            else:# Need to generate the Empty Network here!
                continue
                sub_network = LigandNetwork(edges=set([]), nodes=mols)
                sub_networks.append(sub_network)

        # Connect the Networks:
        con = mst_concatenate(mapper=self.mapper, scorer=self.scorer)
        concat_network = con.concatenate_networks(ligand_networks=sub_networks, n_edges=3)

        return concat_network