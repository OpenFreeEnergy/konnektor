import itertools

from typing import Iterable, Callable, Union, Optional
import functools
from tqdm.auto import tqdm

from gufe import SmallMoleculeComponent, AtomMapper

from openfe.setup.ligand_network import LigandNetwork    # only temproary
from ._abstract_ligand_network_planner import easyLigandNetworkPlanner

from sklearn.cluster import KMeans
from scikit_mol.fingerprints import RDKitFingerprintTransformer, MorganFingerprintTransformer
from scikit_mol.descriptors import MolecularDescriptorTransformer

from konnektor.utils import CompoundDiversityClustering


class DiversityNetworkPlanner(easyLigandNetworkPlanner):
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
        concat_network = planner.concatenate_networks(ligandNetworks=sub_networks, nEdges=3)

        return concat_network