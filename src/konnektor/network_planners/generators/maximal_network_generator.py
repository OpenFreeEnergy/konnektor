# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, Component, LigandNetwork

from .._map_scoring import _score_mappings
from ._abstract_network_generator import NetworkGenerator


class MaximalNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: AtomMapper | Iterable[AtomMapper],
        scorer: Callable[[AtomMapping], float] | None,
        progress: bool = False,
        n_processes: int = 1,
    ):
        r"""
        The ``MaximalNetworkGenerator`` attempts to build a fully connected graph (every node connected to every other node) for given set of ``Component``/s.

        The edges of the graph are ``AtomMapping`` s of pairwise ``Component``\s.
        If not all mappings can be created, it will ignore the mapping failure and return a nearly fully connected graph.

        If multiple ``AtomMapper``/s are provided, but no scorer, the *first valid* ``AtomMapper`` provided will be used.

        .. note::

            This approach is not recommended for Free Energy calculations in application cases, as it is very computationally expensive.
            However, this approach is very important, as many other approaches use the Maximal Network as an initial solution,
            then remove edges to achieve the desired design.

        This class is recommended as an ``initial_edge_lister`` for other network generators.
        The ``MaximalNetworkGenerator`` is parallelized and the number of CPUs can be chosen with the ``n_processes`` argument.

        Parameters
        ----------
        mappers: AtomMapper | list[AtomMapper]
            AtomMapper(s) to use to define the relationship between two ligands. If more than one AtomMapper is provided, all will be tried to find the
            lowest score for each edges.
        scorer: Callable, optional
            Scoring function that takes an AtomMapping and returns a score in [0,1].
        progress: bool, optional
            If True, displays a progress bar, default False.
        n_processes: int
            Number of processes to use for network generation, default 1.
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
        """Create a network with all possible edges.

        Construct a maximally-connected ligand network (n_edges up to $N(N-1)/2$).
        There may be fewer actual mappings than this, because when a mapper cannot
        return a mapping for a given pair, there is simply no suggested mapping for that pair.
        This network is typically used as the starting point for other network
        generators (which then optimize based on the scores) or to debug atom
        mappers (to see which mappings the mapper fails to generate).

        Note that if some mappings are not able to be generated, the resulting graph *may* be disconnected.

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
        TypeError
            If inputs are invalid.
        RuntimeError
            If no mappings were able to be generated.
        """

        components = list(components)

        if self.mappers is None:
            raise TypeError(
                "`mappers` must be an AtomMapper or iterable of AtomMappers to generate a maximal network."
            )

        mappings = _score_mappings(
            possible_edges=list(itertools.combinations(components, 2)),
            scorer=self.scorer,
            mappers=self.mappers,
            n_processes=self.n_processes,
            show_progress=self.progress,
        )

        if len(mappings) == 0:
            raise RuntimeError("Could not generate any mapping!")

        network = LigandNetwork(edges=mappings, nodes=components)
        return network
