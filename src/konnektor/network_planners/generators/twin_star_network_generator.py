# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, Component, LigandNetwork

from konnektor.network_planners._networkx_implementations import RadialNetworkAlgorithm

from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator


class TwinStarNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
        scorer: Callable[[AtomMapping], float] | None,
        n_centers: int = 2,
        n_processes: int = 1,
        progress: bool = False,
        _initial_edge_lister: NetworkGenerator | None = None,
    ):
        """
        The Twin Star Network is an expansion to the Star Network.
        It can be described as multiple star networks that are overlayed.

        The algorithm first calculates all initial edges using `_initial_edge_lister`.
        Next, the in average `n_centers` (default: 2) best performing `Component` s over all
        AtomMapping scores are selected as the central nodes.
        Finally all components are connected to the selected centers, resulting
        in $n_{Transformations} = n_{centers}*(n_{Componentes}-n_{centers})$

        This approach has, in the default version, double the number of edges
        compared to the Star Network, and therefore also an increased graph cost.
        Since the node connectivity is centralized around the `n_centers` nodes,
        the selection of the central ligands is very important, as they have a largeimpact on the overall quality of the network.
        The `n_centers` option allows you to change the Twin Star to a Triplet Star Network or more.

        Parameters
        ----------
        mappers : AtomMapper | list[AtomMapper]
            AtomMapper(s) to use to define the relationship between two ligands. If more than one AtomMapper is provided, all will be tried to find the
            lowest score for each edges.
        scorer : Callable, optional
            Scoring function that takes an AtomMapping and returns a score in [0,1].
        n_centers: int, optional
            Number of central nodes in the network. (default: 2)
        n_processes: int, optional
            Number of processes that to be used for network generation. (default: 1)
        progress: bool, optional
            If True, displays a progress bar, default False.
        _initial_edge_lister:  NetworkGenerator | None
            The NetworkGenerator to use  if the NetworkGenerator requires an initial set of edges.
            If None, a ``MaximalNetworkGenerator`` will be used with the provided mappers and scorer.
        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(
                mappers=mappers, scorer=scorer, n_processes=n_processes
            )

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=RadialNetworkAlgorithm(n_centers=n_centers),
            n_processes=n_processes,
            progress=progress,
            _initial_edge_lister=_initial_edge_lister,
        )

        if isinstance(n_centers, int) and n_centers > 1:
            self.n_centers = n_centers
        else:
            raise ValueError(
                "THe value for n_centers must be an integer > 1. got: " + str(n_centers)
            )

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
        """
        generate a twin star map network for the given compounds.

        Parameters
        ----------
        components: Iterable[Component]
            the components to be used for the LigandNetwork

        Returns
        -------
        LigandNetwork
            a star like network.
        """
        components = list(components)

        # Full Graph Construction
        initial_network = self._initial_edge_lister.generate_ligand_network(components=components)
        mappings = initial_network.edges

        # Translate Mappings to graphable:
        edge_map = {
            (components.index(m.componentA), components.index(m.componentB)): m for m in mappings
        }
        edges = list(sorted(edge_map.keys()))
        weights = [edge_map[k].annotations["score"] for k in edges]

        rg = self.network_generator.generate_network(edges=edges, weights=weights)
        selected_mappings = [edge_map[k] for k in rg.edges]

        return LigandNetwork(edges=selected_mappings, nodes=components)
