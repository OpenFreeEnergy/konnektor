# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import functools
import itertools
from typing import Iterable, Union

from gufe import AtomMapper
from gufe import LigandNetwork, Component
from tqdm.auto import tqdm

from ._abstract_network_generator import NetworkGenerator
from ._parallel_mapping_pattern import _parallel_map_scoring


class MaximalNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: Union[AtomMapper, list[AtomMapper]],
        scorer,
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
            if self.progress is True:
                progress = functools.partial(
                    tqdm, total=total, delay=1.5, desc="Mapping"
                )
            else:
                progress = lambda x: x

            mappings = []
            for component_pair in progress(itertools.combinations(components, 2)):
                best_score = 0.0
                best_mapping = None
                molA = component_pair[0]
                molB = component_pair[1]

                for mapper in self.mappers:
                    try:
                        mapping_generator = mapper.suggest_mappings(molA, molB)
                    except:
                        continue

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
                            print("warning")
                            continue

                if best_mapping is not None:
                    mappings.append(best_mapping)

        if len(mappings) == 0:
            raise RuntimeError("Could not generate any mapping!")

        network = LigandNetwork(edges=mappings, nodes=components)
        return network
