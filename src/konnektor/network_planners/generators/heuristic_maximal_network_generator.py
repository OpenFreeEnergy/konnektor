# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
from collections.abc import Callable, Iterable

import numpy as np
from gufe import AtomMapper, AtomMapping, Component, LigandNetwork

from .._map_scoring import _score_mappings
from ._abstract_network_generator import NetworkGenerator


class HeuristicMaximalNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
        scorer: Callable[[AtomMapping], float] | None,
        n_samples: int = 100,
        progress: bool = False,
        n_processes: int = 1,
    ):
        """The ``HeuristicMaximalNetworkGenerator`` builds, for given set of `Component`/s, a network with `n_samples` edges per `Component`, assuming that each `Component` can be connected to another.

        If a mapping is failed to be generated, the corresponding edge will not be generated, and it is possible that a disconnected graph will be returned.

        This class is recommended for use as an ``initial_edge_lister`` for very large networks.

        .. note ::

            The `HeuristicMaximalNetworkGenerator` can be parallelized and the number of CPUs can be given with  `n_processes`.
            All other approaches in konnektor benefit from this parallelization,ß and you can use this parallelization with `n_processes` key word during class construction.

        Parameters
        ----------
        mappers : AtomMapper | list[AtomMapper]
            AtomMapper(s) to use to define the relationship between two ligands.
        scorer : Callable[[AtomMapping], float] | None
            Scoring function that takes in an atom mapping and returns a score in [0,1].
        n_samples : int, optional
            Max number of of random edges to generate per node, by default 100
        progress : bool, optional
            If True, a progress bar will be displayed, by default False
        n_processes : int, optional
            Number of processes to use for network generation, by default 1
        """

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            n_processes=n_processes,
            network_generator=None,
            progress=progress,
            _initial_edge_lister=self,
        )
        self.n_samples = n_samples

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
        """Create a network with up to n randomly selected edges for each edge.

        Parameters
        ----------
        components : Iterable[Component]
            ``Component``/s to include as nodes in the ``LigandNetwork``.

        Returns
        -------
        LigandNetwork
            ``LigandNetwork`` containing all possible edges, ideally a fully connected graph.

        Raises
        ------
        RuntimeError
            If no mappings were able to be generated.
        """
        components = list(components)

        # Parallel or not Parallel:
        # generate combinations to be searched.
        if len(components) > self.n_samples:
            sample_combinations = []
            for n in components:
                sample_indices = np.random.choice(
                    range(len(components)), size=self.n_samples, replace=False
                )
                sample_combinations.extend(
                    [(n, components[i]) for i in sample_indices if n != components[i]]
                )
        else:
            sample_combinations = list(itertools.combinations(components, 2))

        # todo: what to do if not connected?

        mappings = _score_mappings(
            possible_edges=sample_combinations,
            scorer=self.scorer,
            mappers=self.mappers,
            n_processes=self.n_processes,
            show_progress=self.progress,
        )

        if len(mappings) == 0:
            raise RuntimeError("Could not generate any mapping!")

        network = LigandNetwork(mappings, nodes=components)
        return network
