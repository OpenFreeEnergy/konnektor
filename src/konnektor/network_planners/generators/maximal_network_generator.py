# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import functools
import itertools
from typing import Iterable

from gufe import AtomMapper
from gufe import LigandNetwork, Component
from tqdm.auto import tqdm

from ._abstract_network_generator import NetworkGenerator
from ._parallel_mapping_pattern import _parallel_map_scoring


class MaximalNetworkGenerator(NetworkGenerator):
    def __init__(self, mapper: AtomMapper, scorer, progress: bool = False,
                 n_processes: int = 1):
        """
        The `MaximalNetworkGenerator` builds for given set of `Component`s a fully connected graph under the assumption each `Component` can be connected to another.
        The `Transformation`s of this graph are realized as `AtomMapping`s of pairwise `Component`s. If not all mappings can be created, it will ignore the mapping failure, and return a nearly fully connected graph.

        Note: This approach is not very suitable for Free Energy calculations in application cases. However, this approach is very important, as all above approaches use this as an initial solution, they filter down to gain the desired design.


        This class is recommended as initial_edge_lister for other approaches.
        > **Note**: the `MaximalNetworkGenerator` is parallelized and the number of CPUs can be given with  `n_processes`. 
        > All other approaches in Konnektor benefit from this parallelization and you can use this parallelization with `n_processes` key word during class construction.
        
        Parameters
        ----------
        mapper: AtomMapper
            the atom mapper is required, to define the connection between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        n_processes: int
            number of processes that can be used for the network generation. (default: 1)
        """

        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=None,
                         n_processes=n_processes,
                         _initial_edge_lister=self)
        self.progress = progress

    def generate_ligand_network(self, components: Iterable[
        Component]) -> LigandNetwork:
        """Create a network with all possible proposed mappings.

        This will attempt to create (and optionally score) all possible mappings
        (up to $N(N-1)/2$ for each mapper given). There may be fewer actual
        mappings that this because, when a mapper cannot return a mapping for a
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
        if (self.n_processes > 1):
            mappings = _parallel_map_scoring(
                possible_edges=itertools.combinations(
                    components, 2),
                scorer=self.scorer,
                mapper=self.mapper,
                n_processes=self.n_processes,
                show_progress=self.progress)
        else:  # serial variant
            if self.progress is True:
                progress = functools.partial(tqdm, total=total, delay=1.5,
                                             desc="Mapping")
            else:
                progress = lambda x: x

            mapping_generator = itertools.chain.from_iterable(
                self.mapper.suggest_mappings(molA, molB)
                for molA, molB in
                progress(itertools.combinations(components, 2))
            )
            if self.scorer:
                mappings = [
                    mapping.with_annotations({'score': self.scorer(mapping)})
                    for mapping in mapping_generator]
            else:
                mappings = list(mapping_generator)

        network = LigandNetwork(mappings, nodes=components)
        return network
