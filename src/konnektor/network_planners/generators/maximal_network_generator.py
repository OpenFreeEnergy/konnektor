# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
from collections.abc import Callable, Iterable

from gufe import AtomMapper, Component, LigandNetwork

from .._map_scoring import _parallel_map_scoring, _serial_map_scoring
from ._abstract_network_generator import NetworkGenerator


class MaximalNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
        scorer: Callable | None,
        progress: bool = False,
        n_processes: int = 1,
    ):
        """
        The ``MaximalNetworkGenerator`` attempts to build a fully connected graph (every node connected to every other node) for given set of `Component`/s.

        The edges of the graph are ``Transformation`` s, which contain ``AtomMapping`` s of pairwise ``Component``/s.
        If not all mappings can be created, it will ignore the mapping failure and return a nearly fully connected graph.

        If multiple ``AtomMapper``/s are provided, but no scorer, the *first valid* ``AtomMapper`` provided will be used.

        ... note::
        This approach is not recommended for Free Energy calculations in application cases, as it is very computationally expensive.
        However, this approach is very important, as all other approaches use the Maximal Network as an initial solution,
        then remove edges to achieve the desired design.

        This class is recommended as an ``initial_edge_lister`` for other network generators.
        The ``MaximalNetworkGenerator`` is parallelized and the number of CPUs can be chosen with the ``n_processes`` argument.

        Parameters
        ----------
        mappers: Union[AtomMapper, list[AtomMapper]]
            ``AtomMapper`` to use to define the relationship between two ligands.
        scorer: Callable, optional
            Scoring function that takes in an atom mapping and returns a score in [0,1].
        progress: bool, optional
            If True, a progress bar will be displayed. (default: False)
        n_processes: int
            Number of processes to use for network generation. (default: 1)
        """

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=None,
            n_processes=n_processes,
            progress=progress,
            _initial_edge_lister=self,
        )

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
        """Create a network with all possible proposed mappings.

        This will attempt to create (and optionally score) all possible mappings
        (up to $N(N-1)/2$ for each mapper given). There may be fewer actual
        mappings than this, because when a mapper cannot return a mapping for a
        given pair, there is simply no suggested mapping for that pair.
        This network is typically used as the starting point for other network
        generators (which then optimize based on the scores) or to debug atom
        mappers (to see which mappings the mapper fails to generate).

        Parameters
        ----------
        components : Iterable[SmallMoleculeComponent]
            ``SmallMoleculeComponent``/s to include as nodes in the ``LigandNetwork``.

        Returns
        -------
        LigandNetwork
            ``LigandNetwork`` containing all possible mappings, ideally a fully connected graph.
        """

        components = list(components)
        total = len(components) * (len(components) - 1) // 2

        # Parallel or not Parallel:
        if self.n_processes > 1:
            mappings = _parallel_map_scoring(
                possible_edges=itertools.combinations(components, 2),
                scorer=self.scorer,
                mappers=self.mappers,
                n_processes=self.n_processes,
                show_progress=self.progress,
            )
        else:  # serial variant
            mappings = _serial_map_scoring(
                possible_edges=itertools.combinations(components, 2),
                scorer=self.scorer,
                mappers=self.mappers,
                n_edges_to_score=total,
                show_progress=self.progress,
            )

        if len(mappings) == 0:
            raise RuntimeError("Could not generate any mapping!")

        # TODO: raise an error or warning if disconnected?

        network = LigandNetwork(edges=mappings, nodes=components)
        return network
