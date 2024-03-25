import logging
import inspect
from typing import Iterable, Callable, Union

from tqdm import tqdm

# Clustering
from sklearn.cluster import KMeans
from scikit_mol.fingerprints import RDKitFingerprintTransformer, MorganFingerprintTransformer
from scikit_mol.descriptors import MolecularDescriptorTransformer

from gufe import SmallMoleculeComponent, LigandNetwork
from ...network_tools import cluster_compound, append_node, concatenate_networks
from .. import MinimalSpanningTreeLigandNetworkPlanner, CyclicLigandNetworkPlanner
from ..concatenator import MstConcatenate
from ._abstract_ligand_network_planner import LigandNetworkPlanner

log = logging.getLogger()


class DiversityNetworkPlanner(LigandNetworkPlanner):
    def __init__(self, mapper, scorer,
                 node_featurizer=RDKitFingerprintTransformer(),
                 clustering=KMeans(n_clusters=3),
                 sub_network_planners: Iterable[LigandNetworkPlanner] = (
                         CyclicLigandNetworkPlanner, MinimalSpanningTreeLigandNetworkPlanner),
                 concatenator: MstConcatenate = MstConcatenate,
                 ):
        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=None)
        self.feat = node_featurizer
        self.cluster = clustering

        self.sub_network_planners = [c(mapper=mapper, scorer=scorer) if inspect.isclass(c) else c
                                     for c in sub_network_planners]
        self.concatenator = concatenator(mapper=mapper, scorer=scorer) if inspect.isclass(
            concatenator) else concatenator

    def generate_ligand_network(self,
                                nodes: Iterable[SmallMoleculeComponent],
                                progress: Union[bool, Callable[[Iterable], Iterable]] = True,
                                # allow_disconnected=True
                                ):

        # Step 1: Seperate nodes by diversity
        self.clusters = cluster_compound(compounds=nodes,
                                         featurize=self.feat,
                                         cluster=self.cluster)

        # Step 2:  Sub Network, based on clusters

        self.sub_networks = []
        for cID, mols in tqdm(self.clusters.items(), desc="Build Cluster Networks"):
            if (cID >= 0):  # Noise cluster is not for subnetworks
                if (len(mols) > 1):
                    for network_planner in self.sub_network_planners:
                        try:
                            sub_network = network_planner(mols)
                            self.sub_networks.append(sub_network)
                            break
                        except Exception as err:
                            print("ERR", "\n".join(err.args))
                            continue
                else:  # Need to generate the Empty Network here!
                    sub_network = LigandNetwork(edges=set([]), nodes=mols)
                    self.sub_networks.append(sub_network)

        # step 3: Connect the Networks:
        log.info("Concatenate Networks")
        concat_network = concatenate_networks(networks=self.sub_networks,
                                              concatenator=self.concatenator)

        # step 4: has the clustering a noise cluster
        if -1 in self.clusters:
            for mol in tqdm(self.clusters[-1], desc="add Noise Mols"):
                concat_network = append_node(network=concat_network, compound=mol,
                                             concatenator=self.concatenator)

        return concat_network
