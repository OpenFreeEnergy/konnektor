# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import logging
import warnings
from collections.abc import Iterable

from gufe import AtomMapper, LigandAtomMapping, LigandNetwork

from ._abstract_network_concatenator import NetworkConcatenator
from .mst_concatenator import MstConcatenator

log = logging.getLogger(__name__)


class RedundantMstConcatenator(MstConcatenator):
    def __init__(
        self,
        mappers: AtomMapper | Iterable[AtomMapper] | None,
        scorer,
        n_redundancy: int = 2,
        n_processes: int = 1,
        _initial_edge_lister: NetworkConcatenator | None = None,
    ):
        """
        A NetworkConcatenator that connects subnetworks with `n_redundancy`
        overlaid minimum spanning trees.

        Each pass excludes ligand-pair bridges selected in earlier passes.
        Only complete spanning trees are retained.

        Parameters
        ----------
        mappers: AtomMapper
            AtomMapper(s) to use to propose mappings.
            If more than one AtomMapper is provided, the mapping with the best score (as scored by `scorer`) will be used.
        scorer: Callable[[AtomMapping], float] | None
            Callable which takes a AtomMapping and returns a float in [0,1].
        n_redundancy: int, optional
            Number of spanning trees to overlay, by default 2.
        n_processes: int, optional
            Number of processes that can be used for the network generation, by default 1.
        """
        super().__init__(
            mappers=mappers,
            scorer=scorer,
            n_processes=n_processes,
            _initial_edge_lister=_initial_edge_lister,
        )
        if n_redundancy < 1:
            raise ValueError(f"n_redundancy must be at least 1, got {n_redundancy}")
        self.n_redundancy = n_redundancy

    def _concatenate_networks(
        self,
        ligand_networks: list[LigandNetwork],
        exclude: set[frozenset],
    ) -> LigandNetwork:
        """
        Concatenate the given networks with `n_redundancy` overlaid spanning trees.

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

        Raises
        ------
        RuntimeError
            If any input LigandNetwork is disconnected.
        """
        disconnected_inputs = [n for n in ligand_networks if not n.is_connected()]
        if disconnected_inputs:
            raise RuntimeError(
                f"{len(disconnected_inputs)} of {len(ligand_networks)} input "
                f"subnetworks are disconnected. "
                f"RedundantMstConcatenator expects connected LigandNetworks; "
                f"use connected_subnetworks to split a disconnected network first."
            )

        if len(ligand_networks) == 1:
            return ligand_networks[0]

        # Overlay n_redundancy spanning trees, excluding edges from earlier passes
        bridges: list[LigandAtomMapping] = []
        for n in range(self.n_redundancy):
            new_bridges = self._select_mst_bridges(ligand_networks, exclude)
            if not new_bridges:
                if n > 0:
                    warnings.warn(
                        f"Could only build {n} redundant spanning tree(s) of the "
                        f"{self.n_redundancy} requested; could not form another "
                        f"complete spanning tree from the remaining mappable edges."
                    )
                break
            bridges.extend(new_bridges)
            exclude |= {frozenset((b.componentA, b.componentB)) for b in new_bridges}

        return self._build_concatenated_network(ligand_networks, bridges)
