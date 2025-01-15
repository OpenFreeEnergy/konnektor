# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import functools
import inspect
import logging
from typing import Iterable, Union

from gufe import AtomMapper, Component, LigandNetwork

# Clustering
from scikit_mol.fingerprints import (
    MorganFingerprintTransformer,
    RDKitFingerprintTransformer,
)
from sklearn.cluster import HDBSCAN, KMeans
from tqdm import tqdm

from konnektor.network_tools.clustering.component_diversity_clustering import (
    ComponentsDiversityClusterer,
)

from ...network_tools import append_component, concatenate_networks
from ...network_tools.clustering._abstract_clusterer import _AbstractClusterer
from ..concatenators import MstConcatenator
from ..concatenators._abstract_network_concatenator import NetworkConcatenator
from ._abstract_network_generator import NetworkGenerator
from .cyclic_network_generator import CyclicNetworkGenerator
from .star_network_generator import StarNetworkGenerator

log = logging.getLogger()
log.setLevel(logging.INFO)


class ClusteredNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        sub_network_planners: Iterable[NetworkGenerator] = (CyclicNetworkGenerator,),
        concatenator: NetworkConcatenator = MstConcatenator,
        clusterer: ComponentsDiversityClusterer = ComponentsDiversityClusterer(
            featurize=RDKitFingerprintTransformer(), cluster=KMeans(n_clusters=3)
        ),
        mappers: Union[AtomMapper, list[AtomMapper]] = None,  # include None in this union?
        scorer=None,
        n_processes: int = 1,
        progress: bool = False,
    ):
        """
        Implements the general concept of nd-space clustered networks and provides the logic.

            The algorithm works as follows:
            1. Cluster `Component` s with the `clusterer` obj.
            2. Build sub-networks in the clusters using the `sub_network_planners`.
            3. Concatenate all sub-networks using the `concatenator` to build the final network.


        Parameters
        ----------
        sub_network_planners: Iterable[NetworkGenerator]
            NetworkGenerator(s) used to translate clusters to sub-networks.
        concatenator: NetworkConcatenator
            A `NetworkConcatenator` used to connect sub-networks.
        clusterer: ComponentsDiversityClusterer
            Separates the `Component` s along the first dimension.
        mappers:  Union[AtomMapper, list[AtomMapper]]
            Defines the connection between two ligands if `NetworkConcatenator` s or  `NetworkGenerator` s are provided. Otherwise, (?) (default:None)
        scorer: AtomMappingScorer
            scoring function evaluating an `AtomMapping`, and giving a score between [0,1], if only `NetworkConcatenator` or `NetworkGenerator` classes are passed
        progress: bool, optional
            if True a progress bar will be displayed. (default: False)
        n_processes: int
            number of processes that can be used for the network generation. (default: 1)

        """

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=None,
            progress=progress,
            n_processes=n_processes,
        )
        self.clusterer = clusterer

        if hasattr(self.clusterer.cluster, "n_jobs"):
            self.clusterer.cluster.njobs = n_processes
        if hasattr(self.clusterer.featurize, "n_jobs"):
            self.clusterer.featurize.njobs = n_processes

        if not isinstance(sub_network_planners, (tuple, list)):
            sub_network_planners = [sub_network_planners]

        self.sub_network_planners = []
        for sub_net_planner in sub_network_planners:
            if inspect.isclass(sub_net_planner):
                sub_net_planner_obj = sub_net_planner(mappers=mappers, scorer=scorer)
            else:
                sub_net_planner_obj = sub_net_planner

            sub_net_planner_obj.n_processes = n_processes
            self.sub_network_planners.append(sub_net_planner_obj)

        self.concatenator = (
            concatenator(mappers=mappers, scorer=scorer)
            if inspect.isclass(concatenator)
            else concatenator
        )
        self.concatenator.n_processes = n_processes
        self.progress = progress

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
        """
        Create a network with n randomly selected edges for possible proposed mappings.

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
        log.info("Clusters: " + str(self.clusters))

        # Step 2:  Sub Network, based on clusters
        log.info("Build Sub-Networks")
        self.sub_networks = []
        if self.progress is True:
            progress = functools.partial(
                tqdm, total=len(self.clusters), delay=1.5, desc="Build Cluster Networks"
            )
        else:
            progress = lambda x: x

        if len(self.clusters) == 1 and -1 in self.clusters:
            for network_planner in self.sub_network_planners:
                try:
                    sub_network = network_planner.generate_ligand_network(self.clusters[-1])
                    break
                except Exception as err:
                    print("ERR", "\n".join(err.args))
                    continue
            self.sub_networks.append(sub_network)
        else:
            for cID, mols in progress(self.clusters.items()):
                if cID >= 0:  # Noise cluster is not for subnetworks
                    if len(mols) > 1:
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
            concat_network = concatenate_networks(
                networks=self.sub_networks, concatenator=self.concatenator
            )

        # step 4: has the clustering a noise cluster
        if -1 in self.clusters:
            if self.progress is True:
                progress = functools.partial(
                    tqdm, total=len(self.clusters[-1]), delay=1.5, desc="add Noise Mols"
                )
            else:
                progress = lambda x: x

            for mol in progress(self.clusters[-1]):
                concat_network = append_component(
                    network=concat_network,
                    component=mol,
                    concatenator=self.concatenator,
                )

        return concat_network


class StarrySkyNetworkGenerator(ClusteredNetworkGenerator):
    def __init__(
        self,
        mappers: Union[AtomMapper, list[AtomMapper]],
        scorer,
        clusterer: _AbstractClusterer = ComponentsDiversityClusterer(
            featurize=MorganFingerprintTransformer(),
            cluster=HDBSCAN(metric="jaccard", min_cluster_size=3, alpha=1 / 2048),
        ),
        n_processes: int = 1,
        progress: bool = False,
    ):
        """
        The StarrySkyNetworkGenerator is an advanced network algorithm,
        that clusters the provided `Component` s and builds up a network from this.

        The approach follows the following steps:
        1. Component clustering:
        a. Translate the Molecules into Morgan Fingerprints. (default)
        b. Cluster the Morgan Fingerprints with HDBSCAN. (default)
        2. Build Sub-Star Networks in each Cluster using the `StarNetworkGenerator`.
        3. Concatenate the Sub-Star Networks to the final Starry  Sky Network, with 3 `Transformations` per cluster pair using the `MSTConcatenator`.

        This approach allows in comparison to the Star Network, to build a network containing multiple centers imopoving the graph score.
        Still adding a limited amount of `Transformation` s increasing the computational cost, but not as much `Transformations` as with the Twin Star Network would be generated.
        So the Starry Sky Network is a compromise betwen graph score optimization and number of `Transformations`.

        Parameters
        ----------
        mapper:  Union[AtomMapper, list[AtomMapper]]
            the atom mapper is required, to define the connection between two 'Component's
        scorer: AtomMappingScorer
            scoring function evaluating an `AtomMapping`, and giving a score between [0,1]
        clusterer: ComponentsDiversityClusterer
            This class is seperating the `Component` s along the first dimension.
        progress: bool, optional
            if True a progress bar will be displayed. (default: False)
        n_processes: int
            number of processes that can be used for the network generation. (default: 1)

        """

        super().__init__(
            clusterer=clusterer,
            sub_network_planners=[StarNetworkGenerator],
            concatenator=MstConcatenator,
            mappers=mappers,
            scorer=scorer,
            progress=progress,
            n_processes=n_processes,
        )
