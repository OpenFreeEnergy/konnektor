# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
import logging
from collections.abc import Iterable

from gufe import AtomMapper, LigandAtomMapping, LigandNetwork

from ...network_planners._map_scoring import _score_mappings
from .._networkx_implementations import MstNetworkAlgorithm
from ._abstract_network_concatenator import NetworkConcatenator

log = logging.getLogger(__name__)


class MstConcatenator(NetworkConcatenator):
    def __init__(
        self,
        mappers: AtomMapper | Iterable[AtomMapper] | None,
        scorer,
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
        self, networkA: LigandNetwork, networkB: LigandNetwork
    ) -> list[LigandAtomMapping]:
        """Score every bipartite candidate edge between two subnetworks."""
        possible_edges = [(na, nb) for na in networkA.nodes for nb in networkB.nodes]
        return _score_mappings(
            possible_edges=possible_edges,
            scorer=self.scorer,
            mappers=self.mappers,
            n_processes=self.n_processes,
            show_progress=self.progress,
        )

    def _spanning_tree_pairs(
        self,
        best_mapping_by_pair: dict[tuple[int, int], LigandAtomMapping],
        n_networks: int,
    ) -> list[tuple[int, int]]:
        """
        Build an MST over subnetworks to decide which subnetworks to join.

        Parameters
        ----------
        best_mapping_by_pair : dict[tuple[int, int], LigandAtomMapping]
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
            raise RuntimeError(
                "Could not connect all subnetworks. No mappable edges exist "
                "between the subnetworks."
            )
        subnetwork_edges = list(best_mapping_by_pair)
        # Get the score of the best possible mapping between the subnetworks.
        subnetwork_scores = [
            best_mapping_by_pair[pair].annotations["score"] for pair in subnetwork_edges
        ]
        # Create an MST where each node is a subnetwork
        mst = self.network_generator.generate_network(
            subnetwork_edges, subnetwork_scores, n_edges=n_networks - 1
        )
        if not mst.connected:
            raise RuntimeError(
                "Could not connect all subnetworks. No mappable path exists "
                "between some subnetworks."
            )
        # Reorder the subnetwork indices to match keys in best_mapping_by_pair
        return [(min(i, j), max(i, j)) for i, j in mst.edges]

    def _connect_subnetworks_mst(
        self,
        ligand_networks: list[LigandNetwork],
    ) -> list[LigandAtomMapping]:
        # Find the best scored edge for every pair of subnetworks
        best_mapping_by_pair = {}
        for i, j in itertools.combinations(range(len(ligand_networks)), 2):
            mappings = self._score_pair_edges(ligand_networks[i], ligand_networks[j])
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
        selected_bridges: list[LigandAtomMapping],
    ) -> LigandNetwork:
        # Add the original subnetworks
        edges = list(selected_bridges)
        nodes = set()
        for network in ligand_networks:
            edges.extend(network.edges)
            nodes.update(network.nodes)

        concat_network = LigandNetwork(edges=edges, nodes=nodes)

        return concat_network

    def concatenate_networks(self, ligand_networks: Iterable[LigandNetwork]) -> LigandNetwork:
        """
        Concatenate the given networks.

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            LigandNetworks to concatenate.

        Returns
        -------
        LigandNetwork
            The concatenated LigandNetwork.
        """

        ligand_networks = list(ligand_networks)
        if not ligand_networks:
            raise ValueError("At least one LigandNetwork is required")

        disconnected_inputs = [i for i, n in enumerate(ligand_networks) if not n.is_connected()]
        if disconnected_inputs:
            raise RuntimeError(
                f"Input subnetworks {disconnected_inputs} are disconnected. "
                f"MstConcatenator expects connected LigandNetworks; "
                f"use decompose_network to split a disconnected network first."
            )

        log.info(
            f"Number of edges in individual networks:\n"
            f"{sum(len(s.edges) for s in ligand_networks)}/"
            f"{[len(s.edges) for s in ligand_networks]}"
        )

        if len(ligand_networks) == 1:
            return ligand_networks[0]

        selected_bridges = self._connect_subnetworks_mst(ligand_networks)
        concat_network = self._build_concatenated_network(ligand_networks, selected_bridges)

        return concat_network
