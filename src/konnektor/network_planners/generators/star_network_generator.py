# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import functools
from collections.abc import Iterable

from gufe import AtomMapper, Component, LigandNetwork
from tqdm import tqdm
import warnings

from konnektor.network_planners._networkx_implementations import RadialNetworkAlgorithm

from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator


class StarNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
        scorer,
        n_processes: int = 1,
        progress: bool = False,
        _initial_edge_lister: NetworkGenerator = None,
    ):
        """
        The Star Network is one of the most edge efficient layouts, it basically places
        all `Transformations` around one central `Component`.

        The algorithm constructs in a first step all possible `Transformations`.
        Next it selects in the default variant the in average best transformation score performing `Component` as the central component.
        Finally all Components are connected with a `Transformation` to the central `Component`

        The Star Network is most edge efficient, but not most graph score efficient, as it has to find a
        central `Component`, which usually is a compromise for all 'Component's.
        From a robustness point of view, the Star Network, will immediately be disconnected if one `Transformation` fails.
        However the loss of `Component` s is very limited, as only one ligand is lost per `Transformation` failure.

        Parameters
        ----------
        mapper : AtomMapper
            the atom mapper is required, to define the connection between two ligands.
        scorer : AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        n_processes: int, optional
            number of processes that can be used for the network generation. (default: 1)
        progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        _initial_edge_lister: NetworkPlanner, optional
            this NetworkPlanner is used to give the initial set of edges. For standard usage, the Maximal NetworPlanner is used.
            However in large scale approaches, it might be interesting to use the heuristicMaximalNetworkPlanner..
            (default: MaximalNetworkPlanner)
        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(
                mappers=mappers, scorer=scorer, n_processes=n_processes
            )

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=RadialNetworkAlgorithm(),
            n_processes=n_processes,
            progress=progress,
            _initial_edge_lister=_initial_edge_lister,
        )

    def generate_ligand_network(
        self, components: Iterable[Component], central_component: Component = None
    ) -> LigandNetwork:
        """
        generate a star map network for the given compounds. if a central component is defined,
        the planning stage is shortcutted to only connect the ligands to the central component.

        Parameters
        ----------
        components: Iterable[Component]
            the components to be used for the LigandNetwork

        central_component: Component, optional
            the central component can be given, in order to shortcut the calculations and enforce the central ligand.

        Returns
        -------
        LigandNetwork
            a star like network.
        """
        components = list(components)

        # If central Ligand not defined, get it from full graph.
        if central_component is None:
            # Full Graph Construction
            initial_network = self._initial_edge_lister.generate_ligand_network(
                components=components
            )
            mappings = initial_network.edges

            # Translate Mappings to graphable:
            edge_map = {
                (components.index(m.componentA), components.index(m.componentB)): m
                for m in mappings
            }
            edges = list(sorted(edge_map.keys()))
            weights = [edge_map[k].annotations["score"] for k in edges]

            rg = self.network_generator.generate_network(edges=edges, weights=weights)
            selected_mappings = [edge_map[k] for k in rg.edges]

        else:  # Given central ligands: less effort. - Trivial Case
            # TODO: commenting this out because scorer isn't used - should it be used?
            # if self.scorer is None:
            #     scorer = lambda x: -1
            # else:
            #     scorer = self.scorer

            if self.progress is True:
                progress = functools.partial(tqdm, total=len(components), delay=1.5, desc="Mapping")
            else:
                progress = lambda x: x

            selected_mappings = []
            for component in progress(components):
                best_score = 0.0
                best_mapping = None
                molA = central_component
                molB = component

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
                            # TODO: output which mapper is first?
                            warnings.warn("Multiple mappers were provided, but no scorer. Only the first mapper provided will be used.")
                            best_mapping = next(mapping_generator)
                            break
                        except:
                            continue
                if best_mapping is not None:
                    selected_mappings.append(best_mapping)

        return LigandNetwork(edges=selected_mappings, nodes=components)


RadialLigandNetworkPlanner = StarNetworkGenerator
