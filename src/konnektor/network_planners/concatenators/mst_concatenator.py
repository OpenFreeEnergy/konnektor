# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
import logging
from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, LigandNetwork

from .._networkx_implementations import MstNetworkAlgorithm
from ._abstract_network_concatenator import NetworkConcatenator

log = logging.getLogger(__name__)


class MstConcatenator(NetworkConcatenator):
    def __init__(
        self,
        mappers: AtomMapper | Iterable[AtomMapper] | None,
        scorer: Callable[[AtomMapping], float],
        n_processes: int = 1,
        _initial_edge_lister: NetworkConcatenator | None = None,  # TODO: remove this
    ):
        """
        A NetworkConcatenator that connects subnetworks using a minimum
        spanning tree over the subnetworks.

        For each pair of subnetworks, the best-scoring ligand mapping is used
        as the candidate edge between those subnetworks.

        Parameters
        ----------
        mappers: AtomMapper
            AtomMapper(s) to use to propose mappings.
            If more than one AtomMapper is provided, the mapping with the best score (as scored by `scorer`) will be used.
        scorer: Callable[[AtomMapping], float]
            Callable which takes a AtomMapping and returns a float in [0,1].
        n_processes: int, optional
            Number of processes that can be used for the network generation, by default 1.
        """
        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=MstNetworkAlgorithm(),
            n_processes=n_processes,
            _initial_edge_lister=_initial_edge_lister,
        )

    def _select_spanning_tree_pairs(
        self,
        best_mapping_by_network_pair: dict[tuple[int, int], AtomMapping],
        n_networks: int,
    ) -> list[tuple[int, int]]:
        """
        Select the subnetwork pairs that form a complete spanning tree.

        Parameters
        ----------
        best_mapping_by_network_pair: dict[tuple[int, int], AtomMapping]
            The best-scoring mapping for each pair of subnetworks.
            Each key (i, j) identifies ligand_networks[i] and ligand_networks[j].
        n_networks: int
            Number of input subnetworks.

        Returns
        -------
        list[tuple[int, int]]
            Pairs of input subnetwork indices defining the edges of the MST.
        """
        if not best_mapping_by_pair:
            return []
        subnetwork_edges = list(best_mapping_by_network_pair)
        # Get the score of the best possible mapping between the subnetworks.
        subnetwork_scores = [
            best_mapping_by_network_pair[pair].annotations["score"] for pair in subnetwork_edges
        ]
        # Create an MST where each node is a subnetwork
        mst = self.network_generator.generate_network(
            subnetwork_edges,
            subnetwork_scores,
        )

        # Reorder the subnetwork indices to match keys in best_mapping_by_pair
        subnetwork_pairs = [(min(i, j), max(i, j)) for i, j in mst.edges]

        # A spanning tree over n subnetworks must contain exactly n - 1 edges.
        # Discard incomplete spanning forests.
        if len(subnetwork_pairs) != n_networks - 1:
            return []

        return subnetwork_pairs

    def _select_mst_bridges(
        self,
        ligand_networks: list[LigandNetwork],
        exclude: set[frozenset],
    ) -> list[AtomMapping]:
        # Find the best scored edge for every pair of subnetworks
        best_mapping_by_network_pair = {}
        for i, j in itertools.combinations(range(len(ligand_networks)), 2):
            mappings = self._score_inter_network_edges(
                ligand_networks[i],
                ligand_networks[j],
                exclude,
            )
            if mappings:
                # Best ligand mapping for this pair of subnetworks
                best_mapping_by_network_pair[(i, j)] = max(
                    mappings,
                    key=lambda mapping: (
                        mapping.annotations["score"],
                        mapping.key,  # deterministic tie-break
                    ),
                )

        # Identify which subnetworks to connect (MST over subnetworks)
        subnetwork_pairs = self._select_spanning_tree_pairs(best_mapping_by_network_pair, len(ligand_networks))

        # Connect each subnetwork pair with best scored edge
        selected_bridges = [best_mapping_by_network_pair[pair] for pair in subnetwork_pairs]
        return selected_bridges

    def _concatenate_networks(
        self,
        ligand_networks: list[LigandNetwork],
        exclude: set[frozenset],
    ) -> LigandNetwork:
        """
        Concatenate the given ligand_networks, treating each ligand_network as a node in a Minimal Spanning Tree.

        Parameters
        ----------
        ligand_networks: list[LigandNetwork]
            LigandNetworks to concatenate.
        exclude : set[frozenset]
            Unordered ligand pairs that must not be proposed as new connections.

        Returns
        -------
        LigandNetwork
            The concatenated LigandNetwork.
        """
        selected_bridges = self._select_mst_bridges(ligand_networks, exclude)

        return self._assemble_concatenated_network(ligand_networks, selected_bridges)
