import functools
import inspect
import logging
from tqdm import tqdm
from typing import Iterable

from gufe import Component, LigandNetwork, AtomMapper, AtomMappingScorer

# Clustering
from scikit_mol.fingerprints import RDKitFingerprintTransformer, MorganFingerprintTransformer
from sklearn.cluster import HDBSCAN

from ._abstract_network_generator import NetworkGenerator

from .. import RadialLigandNetworkPlanner
from ..concatenator import MstConcatenate
from ...network_tools import append_node, concatenate_networks
from ...network_tools.cluster_components import ComponentsDiversityClustering

log = logging.getLogger()
log.setLevel(logging.INFO)

# Todo: go over this again.

class TwoDimensionalNetworkGenerator(NetworkGenerator):
    def __init__(self,
                 sub_network_planners: Iterable[NetworkGenerator] = (RadialLigandNetworkPlanner),
                 concatenator: MstConcatenate = MstConcatenate,
                 clusterer: ComponentsDiversityClustering = ComponentsDiversityClustering(
                     featurize=RDKitFingerprintTransformer(), cluster=HDBSCAN()),
                 mapper: AtomMapper = None, scorer: AtomMappingScorer = None,
                 nprocesses: int = 1, progress: bool = False
                 ):
        ''' Implements the general concept of multidimensional networks.

        Parameters
        ----------
        clusterer: ComponentsDiversityClustering
            This class is seperating the Components along the first dimension.
        sub_network_planners: Iterable[NetworkGenerator]
            The clusters, are then seperatley translated to sub networks by the sub_network_planners
        concatenator: MstConcatenate
            The concatenator is connecting the different sub networks.
        mapper: AtomMapper
            the atom mapper is required, to define the connection between two ligands, if only concatenator or ligandPlanner classes are passed
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1], if only concatenator or ligandPlanner classes are passed
        progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        nprocesses: int
            number of processes that can be used for the network generation. (default: 1)

        '''

        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=None)
        self.clusterer = clusterer

        if hasattr(self.clusterer.cluster, "n_jobs"):
            self.clusterer.cluster.njobs = nprocesses

        self.sub_network_planners = []
        for sub_net_planner in sub_network_planners:
            if inspect.isclass(sub_net_planner):
                sub_net_planner_obj = sub_net_planner(mapper=mapper, scorer=scorer)
            else:
                sub_net_planner_obj = sub_net_planner
            sub_net_planner_obj.nprocesses = nprocesses
            self.sub_network_planners.append(sub_net_planner_obj)

        self.concatenator = concatenator(mapper=mapper, scorer=scorer) if inspect.isclass(
            concatenator) else concatenator
        self.concatenator.nprocesses = nprocesses
        self.progress = progress

    def generate_ligand_network(self,
                                components: Iterable[Component]
                                ) -> LigandNetwork:
        """Create a network with n randomly selected edges for possible proposed mappings.

        Parameters
        ----------
        components : Iterable[Component]
          the ligands to include in the LigandNetwork

        Returns
        -------
        LigandNetwork
            a complex network.
        """

        # Step 1: Seperate nodes by diversity
        log.info("Clustering")
        self.clusters = self.clusterer.cluster_compounds(components)
        log.info("Clusters: "+str(self.clusters))

        if len(self.clusters) == 1 and -1 in self.clusters:
            print(self.clusters.keys())
            raise ValueError("Found only noise cluster nodes!")

        # Step 2:  Sub Network, based on clusters
        log.info("Build Sub-Networks")
        self.sub_networks = []
        if self.progress is True:
            progress = functools.partial(tqdm, total=len(self.clusters), delay=1.5,
                                         desc="Build Cluster Networks")
        else:
            progress = lambda x: x

        for cID, mols in progress(self.clusters.items()):
            if (cID >= 0):  # Noise cluster is not for subnetworks
                if (len(mols) > 1):
                    for network_planner in self.sub_network_planners:
                        try:
                            sub_network = network_planner.generate_ligand_network(mols)
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
        if len(self.sub_networks) == 1:
            concat_network = self.sub_networks[0]
        else:
            concat_network = concatenate_networks(networks=self.sub_networks,
                                              concatenator=self.concatenator)

        # step 4: has the clustering a noise cluster
        if -1 in self.clusters:
            if self.progress is True:
                progress = functools.partial(tqdm, total=len(self.clusters[-1]), delay=1.5,
                                             desc="add Noise Mols")
            else:
                progress = lambda x: x

            for mol in progress(self.clusters[-1]):
                concat_network = append_node(network=concat_network, compound=mol,
                                             concatenator=self.concatenator)

        return concat_network

class StarrySkyNetworkGenerator(TwoDimensionalNetworkGenerator):
    def __init__(self,
                 clusterer: ComponentsDiversityClustering = ComponentsDiversityClustering(
                     featurize=MorganFingerprintTransformer(), cluster=HDBSCAN(metric="jaccard", alpha=2048)),
                 mapper: AtomMapper = None,
                 scorer: AtomMappingScorer = None,
                 nprocesses: int = 1, progress: bool = False
                 ):
        '''  The StarrySkyNetworkGenerator is an advanced network algorithm,
        that clusters the provided ligands, based on the clusterer class.
        Next the clusters are centered around one ligand into a Star map.
        The star maps are connected with a MSTconcatenate.

        Parameters
        ----------
        clusterer: ComponentsDiversityClustering
            This class is seperating the Components along the first dimension.
        mapper: AtomMapper
            the atom mapper is required, to define the connection between two ligands, if only concatenator or ligandPlanner classes are passed
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1], if only concatenator or ligandPlanner classes are passed
        progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        nprocesses: int
            number of processes that can be used for the network generation. (default: 1)

        '''

        super().__init__(clusterer=clusterer,
                         sub_network_planners=[RadialLigandNetworkPlanner],
                         concatenator = MstConcatenate,
                         mapper=mapper, scorer=scorer,
                         progress=progress,
                         nprocesses=nprocesses,
                         )