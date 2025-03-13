# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, Component, LigandNetwork

from .._map_scoring import _parallel_map_scoring, _serial_map_scoring
from ._abstract_network_generator import NetworkGenerator


class MaximalNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
        scorer: Callable[[AtomMapping], float] | None,
        progress: bool = False,
        n_processes: int = 1,
    ):
        """
        The `MaximalNetworkGenerator` builds a fully connected graph for given set of `Component` s.
        It assumes each `Component` can be connected to every other `Component`.
        The `Transformation` s of this graph are realized as `AtomMapping` s of pairwise `Component` s.
        If not all mappings can be created, it will ignore the mapping failure and return a nearly fully connected graph.

        ... note::
        This approach is not recommended for Free Energy calculations in application cases, as it is very computationally expensive.
        However, this approach is very important, as all other approaches use the Maximal Network as an initial solution,
        then remove edges to achieve the desired design.

        This class is recommended as an `initial_edge_lister` for other approaches.
        The `MaximalNetworkGenerator` is parallelized and the number of CPUs can be given with `n_processes`.
        All other approaches in Konnektor benefit from this parallelization and you can use this parallelization with `n_processes` key word during class construction.

        Parameters
        ----------
        mappers: Union[AtomMapper, list[AtomMapper]]
            the atom mapper is required, to define the connection between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        n_processes: int
            number of processes that can be used for the network generation. (default: 1)
        """

        if isinstance(mappers, list) and scorer is None:
            if len(mappers) > 1:
                raise ValueError("You must provide a scorer when passing in multiple mappers.")

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
          the ligands to include in the LigandNetwork

        Returns
        -------
        LigandNetwork
            a ligand network containing all possible mappings, ideally a fully connected graph.
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

        network = LigandNetwork(edges=mappings, nodes=components)
        return network
