# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import functools
import itertools
from typing import Iterable, Union

import numpy as np
from gufe import AtomMapper, Component, LigandNetwork
from tqdm.auto import tqdm

from ._abstract_network_generator import NetworkGenerator
from ._parallel_mapping_pattern import _parallel_map_scoring

# Todo: is graph connectivity ensured?


class HeuristicMaximalNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: Union[AtomMapper, list[AtomMapper]],
        scorer,
        n_samples: int = 100,
        progress: bool = False,
        n_processes: int = 1,
    ):
        """
        The `HeuristicMaximalNetworkGenerator` builds for given set of `Component` s a set of `n_samples` `Transformation` s per `Component` build network under the assumption each `Component` can be connected to another.
        The `Transformations` of this network are realized as `AtomMapping` s of pairwise `Component` s. If not all mappings can be created, it will ignore the mapping failure, and return a nearly fully connected graph.

        This class is can be used as initial_edge_lister, if there is a large set of `Component` s (check network connectivity!)

        This class is recommended as initial_edge_lister for other approaches.
        > **Note**: the `HeuristicMaximalNetworkGenerator` is parallelized and the number of CPUs can be given with  `n_processes`.
        > All other approaches in Konnektor benefit from this parallelization and you can use this parallelization with `n_processes` key word during class construction.


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
        n_processes: int
            number of processes that can be used for the network generation. (default: 1)


        """
        super().__init__(
            mappers=mappers,
            scorer=scorer,
            n_processes=n_processes,
            network_generator=None,
            progress=progress,
            _initial_edge_lister=self,
        )

        self.progress = progress
        self.n_samples = n_samples

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
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
        total = len(components) * self.n_samples

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
            sample_combinations = itertools.combinations(components, 2)

        # todo: what to do if not connected?

        if self.n_processes > 1:
            mappings = _parallel_map_scoring(
                possible_edges=sample_combinations,
                scorer=self.scorer,
                mappers=self.mappers,
                n_processes=self.n_processes,
                show_progress=self.progress,
            )
        else:  # serial variant
            if self.progress is True:
                progress = functools.partial(tqdm, total=total, delay=1.5, desc="Mapping")
            else:
                progress = lambda x: x

            mappings = []
            for component_pair in progress(sample_combinations):
                best_score = 0.0
                best_mapping = None
                molA = component_pair[0]
                molB = component_pair[1]

                for mapper in self.mappers:
                    mapping_generator = mapper.suggest_mappings(molA, molB)

                    if self.scorer:
                        tmp_mappings = [
                            mapping.with_annotations({"score": self.scorer(mapping)})
                            for mapping in mapping_generator
                        ]

                        if len(tmp_mappings) > 0:
                            tmp_best_mapping = min(
                                tmp_mappings, key=lambda m: m.annotations["score"]
                            )

                            if (
                                tmp_best_mapping.annotations["score"] < best_score
                                or best_mapping is None
                            ):
                                best_score = tmp_best_mapping.annotations["score"]
                                best_mapping = tmp_best_mapping
                    else:
                        try:
                            best_mapping = next(mapping_generator)
                        except:
                            continue
                if best_mapping is not None:
                    mappings.append(best_mapping)

        if len(mappings) == 0:
            raise RuntimeError("Could not generate any mapping!")

        network = LigandNetwork(mappings, nodes=components)
        return network
