import functools
import itertools
from typing import Iterable

import numpy as np
from gufe import Component, LigandNetwork, AtomMapper
from tqdm.auto import tqdm

from ._abstract_network_generator import NetworkGenerator
from ._parallel_mapping_pattern import _parallel_map_scoring


# Todo: is graph connectivity ensured?

class HeuristicMaximalNetworkGenerator(NetworkGenerator):
    def __init__(self, mapper: AtomMapper, scorer, n_samples: int = 100,
                 progress: bool = False,
                 nprocesses: int = 1):
        """
        The Heuristic Maximal Network planner builds for given set of compounds a set of edges per node build graph under the assumption each component can be connected to another.
        The edges of this graph are realized as atom mappings of pairwise components. If not all mappings can be created, it will ignore the mapping failure, and return a nearly fully connected graph.

        This class is can be used as initial_edge_lister

        Parameters
        ----------
        mapper: AtomMapper
            the atom mapper is required, to define the connection between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        n_samples: int
            number of random edges per node.
        progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        nprocesses: int
            number of processes that can be used for the network generation. (default: 1)


        """
        super().__init__(mapper=mapper, scorer=scorer,
                         nprocesses=nprocesses,
                         network_generator=None,
                         _initial_edge_lister=self)

        self.progress = progress
        self.n_samples = n_samples

    def generate_ligand_network(self, components: Iterable[
        Component]) -> LigandNetwork:
        """Create a network with n randomly selected edges for possible proposed mappings.

        Parameters
        ----------
        components : Iterable[Component]
          the ligands to include in the LigandNetwork

        Returns
        -------
        LigandNetwork
            a heuristic max network.
        """
        components = list(components)
        total = len(components) * (len(components) - 1) // 2

        # Parallel or not Parallel:
        # generate combinations to be searched.
        if len(components) > self.n_samples:
            sample_combinations = []
            for n in components:
                sample_indices = np.random.choice(range(len(components)),
                                                  size=self.n_samples,
                                                  replace=False)
                sample_combinations.extend(
                    [(n, components[i]) for i in sample_indices if
                     n != components[i]])
        else:
            sample_combinations = itertools.combinations(components, 2)

        # todo: what to do if not connected?

        if (self.nprocesses > 1):
            mappings = _parallel_map_scoring(
                possible_edges=sample_combinations,
                scorer=self.scorer,
                mapper=self.mapper,
                n_processes=self.nprocesses,
                show_progress=self.progress)
        else:  # serial variant
            if self.progress is True:
                progress = functools.partial(tqdm, total=total, delay=1.5,
                                             desc="Mapping")
            else:
                progress = lambda x: x

            mapping_generator = itertools.chain.from_iterable(
                self.mapper.suggest_mappings(molA, molB)
                for molA, molB in progress(sample_combinations)
            )
            if self.scorer:
                mappings = [
                    mapping.with_annotations({'score': self.scorer(mapping)})
                    for mapping in mapping_generator]
            else:
                mappings = list(mapping_generator)

        network = LigandNetwork(mappings, nodes=components)
        return network
