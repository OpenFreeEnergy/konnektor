# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from typing import Iterable, Union

from gufe import Component, LigandNetwork, AtomMapper

from konnektor.network_planners._networkx_implementations import RadialNetworkAlgorithm
from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator


class TwinStarNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: Union[AtomMapper, list[AtomMapper]],
        scorer,
        n_centers: int = 2,
        n_processes: int = 1,
        progress: bool = False,
        _initial_edge_lister: NetworkGenerator = None,
    ):
        """
        The Twin Star Network is an expansion to the Star Network. It can be described as multiple star networks that are overlayed.

        The algorithm first calculates all possible `Transformation` s for all `Component` s.
        Next, the in average `n_centers` (default: 2) best performing `Component` s over all
        transfromation scores are selected and placed into the center of the network.
        Finally all components are connected to the selected centers, resulting
        in $n_{Transformations} = n_{centers}*(n_{Componentes}-n_{centers})$

        This approach has, in the default version, double the number of `Transformations`
        compared to the Star Network, and therefore also an increased graph cost.
        On the plus side, this approach builds many graph cycles, which could be used to estimate the uncertainty of FE calculations.
        Another important aspect is that the node connectivity is centralized around the `n_centers`
        This means that the selection of the central ligands is very important, as they have a large
        impact on the `Transformations`.
        The `n_centers` option allows you to change the Twin Star to a Triplet Star Network or more.

        Parameters
        ----------
        mapper :  Union[AtomMapper, list[AtomMapper]]
            the atom mapper is required, to define the connection between two ligands.
        scorer : AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        n_centers: int, optional
            the number of centers in the network. (default: 2)
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

        return LigandNetwork(edges=selected_mappings, nodes=components)
