# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from tqdm import tqdm
import functools
import itertools
import logging
from typing import Iterable, Union

from gufe import AtomMapper, LigandNetwork

from ._abstract_network_concatenator import NetworkConcatenator
from ..generators._parallel_mapping_pattern import _parallel_map_scoring

log = logging.getLogger(__name__)


class MaxConcatenator(NetworkConcatenator):
    def __init__(
        self,
        mappers: Union[AtomMapper, list[AtomMapper]],
        scorer,
        n_processes: int = 1,
        show_progress: bool = False,
    ):
        """
        This concatenators is connnecting two Networks with all possible
         mappings. This is usually most useful for initial edge listing.

        Parameters
        ----------
        mapper: AtomMapper
            the atom mapper is required, to define the connection
             between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a
            score between [0,1].
        n_connecting_edges: int, optional
            number of connecting edges. (default: 3)
        n_processes: int
            number of processes that can be used for the network generation.
            (default: 1)
        show_progress: bool
            show progress bar
        """

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=None,
            n_processes=n_processes,
        )
        self.progress = show_progress

    def concatenate_networks(
        self, ligand_networks: Iterable[LigandNetwork]
    ) -> LigandNetwork:
        """

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            an iterable of ligand networks, that shall be connected.

        Returns
        -------
        LigandNetwork
            returns a concatenated LigandNetwork object, containing all
             networks and all possible edges, connecting them.

        """

        log.info(
            f"Number of edges in individual networks:\n"
            f"{sum([len(s.edges) for s in ligand_networks])}/"
            f"{[len(s.edges) for s in ligand_networks]}"
        )

        selected_edges = []
        selected_nodes = []
        for ligandNetworkA, ligandNetworkB in itertools.combinations(
            ligand_networks, 2
        ):
            # Generate Full Bipartite Graph
            nodesA = ligandNetworkA.nodes
            nodesB = ligandNetworkB.nodes
            pedges = [(na, nb) for na in nodesA for nb in nodesB]

            if self.n_processes > 1:
                bipartite_graph_mappings = _parallel_map_scoring(
                    possible_edges=pedges,
                    scorer=self.scorer,
                    mappers=self.mappers,
                    n_processes=self.n_processes,
                    show_progress=self.progress,
                )

            else:  # serial variant
                if self.progress is True:
                    progress = functools.partial(
                        tqdm, total=len(pedges), delay=1.5, desc="Mapping Subnets"
                    )
                else:
                    progress = lambda x: x

                bipartite_graph_mappings = []
                for component_pair in progress(pedges):
                    best_score = 0.0
                    best_mapping = None
                    molA = component_pair[0]
                    molB = component_pair[1]

                    for mapper in self.mappers:
                        try:
                            mapping_generator = mapper.suggest_mappings(molA, molB)
                        except:
                            continue

                        if self.scorer:
                            tmp_mappings = [
                                mapping.with_annotations(
                                    {"score": self.scorer(mapping)}
                                )
                                for mapping in mapping_generator
                            ]

                            if len(tmp_mappings) > 0:
                                tmp_best_mapping = min(
                                    tmp_mappings, key=lambda m: m.annotations["score"]
                                )

                                if (
                                    tmp_best_mapping.annotations["score"] < best_score
                                    or best_mapping is None
                                ):
                                    best_score = tmp_best_mapping.annotations["score"]
                                    best_mapping = tmp_best_mapping
                        else:
                            try:
                                best_mapping = next(mapping_generator)
                            except:
                                print("warning")
                                continue

                    if best_mapping is not None:
                        bipartite_graph_mappings.append(best_mapping)

            # Add network connecting edges
            selected_edges.extend(bipartite_graph_mappings)

        # Constructed final Edges:
        # Add all old network edges:
        for network in ligand_networks:
            selected_edges.extend(network.edges)
            selected_nodes.extend(network.nodes)

        concat_LigandNetwork = LigandNetwork(
            edges=selected_edges, nodes=set(selected_nodes)
        )

        log.info(f"Total Concatenated Edges: {len(selected_edges)} ")

        if not concat_LigandNetwork.is_connected():
            raise RuntimeError("could not build a connected network!")

        return concat_LigandNetwork
