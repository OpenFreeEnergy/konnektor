# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
import logging
from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, LigandNetwork

from ...network_planners._map_scoring import _score_mappings
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
        A NetworkConcatenator that connects subnetworks by treating each
        subnetwork as a node in an MST.

        Parameters
        ----------
        mappers: AtomMapper
            AtomMapper(s) to use to propose mappings.
            If more than one AtomMapper is provided, the mapping with the best score (as scored by `scorer`) will be used.
        scorer: Callable[[AtomMapping], float] | None
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

    def _score_pair_edges(
        self,
        networkA: LigandNetwork,
        networkB: LigandNetwork,
        exclude: set[frozenset],
    ) -> list[AtomMapping]:
        """Score every bipartite candidate edge between two subnetworks."""
        possible_edges = self._generate_bipartite_edges(networkA, networkB, exclude)
        return _score_mappings(
            possible_edges=possible_edges,
            scorer=self.scorer,
            mappers=self.mappers,
            n_processes=self.n_processes,
            show_progress=self.progress,
        )

    def _spanning_tree_pairs(
        self,
        best_mapping_by_pair: dict[tuple[int, int], AtomMapping],
        n_networks: int,
    ) -> list[tuple[int, int]]:
        """
        Build an MST over subnetworks to decide which subnetworks to join.

        Parameters
        ----------
        best_mapping_by_pair : dict[tuple[int, int], AtomMapping]
            The best-scoring mapping for each pair of subnetworks.
            Each key (i, j) identifies ligand_networks[i] and ligand_networks[j].
        n_networks : int
            Number of input subnetworks.

        Returns
        -------
        list[tuple[int, int]]
            Pairs of input subnetwork indices defining the edges of the MST.
        """
        if not best_mapping_by_pair:
            return []
        subnetwork_edges = list(best_mapping_by_pair)
        # Get the score of the best possible mapping between the subnetworks.
        subnetwork_scores = [
            best_mapping_by_pair[pair].annotations["score"] for pair in subnetwork_edges
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

    def _connect_subnetworks_mst(
        self,
        ligand_networks: list[LigandNetwork],
        exclude: set[frozenset],
    ) -> list[AtomMapping]:
        # Find the best scored edge for every pair of subnetworks
        best_mapping_by_pair = {}
        for i, j in itertools.combinations(range(len(ligand_networks)), 2):
            mappings = self._score_pair_edges(ligand_networks[i], ligand_networks[j], exclude)
            if mappings:
                best_mapping_by_pair[(i, j)] = max(
                    mappings,
                    key=lambda mapping: (
                        mapping.annotations["score"],
                        mapping.key,  # deterministic tie-break
                    ),
                )

        # Identify which subnetworks to connect (MST over subnetworks)
        subnetwork_pairs = self._spanning_tree_pairs(best_mapping_by_pair, len(ligand_networks))

        # Connect each subnetwork pair with best scored edge
        selected_bridges = [best_mapping_by_pair[pair] for pair in subnetwork_pairs]
        return selected_bridges

    def _build_concatenated_network(
        self,
        ligand_networks: list[LigandNetwork],
        selected_bridges: list[AtomMapping],
    ) -> LigandNetwork:
        # Add the original subnetworks
        edges = list(selected_bridges)
        nodes = set()
        for network in ligand_networks:
            edges.extend(network.edges)
            nodes.update(network.nodes)

        concat_network = LigandNetwork(edges=edges, nodes=nodes)
        if not concat_network.is_connected():
            raise RuntimeError(
                "Could not connect all subnetworks. No mappable edges exist "
                "between some subnetworks (possibly all excluded via avoid_edges)."
            )

        return concat_network

    def concatenate_networks(
        self,
        ligand_networks: Iterable[LigandNetwork],
        exclude_edges: Iterable[AtomMapping] | None = None,
    ) -> LigandNetwork:
        """
        Concatenate the given networks.

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            LigandNetworks to concatenate.
        exclude_edges: Iterable[AtomMapping], optional
            Mappings that cannot be proposed as new connections which is useful for excluding edges
            that had already failed. If excluding these edges leaves the network unbridgeable, an error is raised.
            Default: None

        Returns
        -------
        LigandNetwork
            The concatenated LigandNetwork.

        Raises
        ------
        RuntimeError
            If the network cannot be connected, either because the input network
            was disconnected or no mappable edges could be found between the subnetworks.
        """

        ligand_networks = list(ligand_networks)
        if not ligand_networks:
            raise ValueError("At least one LigandNetwork is required")

        exclude = self._normalize_excluded_edges(exclude_edges)

        disconnected_inputs = [n for n in ligand_networks if not n.is_connected()]
        if disconnected_inputs:
            raise RuntimeError(
                f"{len(disconnected_inputs)} of {len(ligand_networks)} input "
                f"subnetworks are disconnected. "
                f"MstConcatenator expects connected LigandNetworks; "
                f"use connected_subnetworks to split a disconnected network first."
            )

        log.info(
            f"Number of edges in individual networks:\n"
            f"{sum(len(s.edges) for s in ligand_networks)}/"
            f"{[len(s.edges) for s in ligand_networks]}"
        )

        if len(ligand_networks) == 1:
            return ligand_networks[0]

        selected_bridges = self._connect_subnetworks_mst(ligand_networks, exclude)
        concat_network = self._build_concatenated_network(ligand_networks, selected_bridges)

        return concat_network
