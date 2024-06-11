# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from typing import Iterable

from gufe import Component, LigandNetwork, AtomMapper

from konnektor.network_planners._networkx_implementations import \
    RadialNetworkAlgorithm
from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator


class TwinStarNetworkGenerator(NetworkGenerator):

    def __init__(self, mapper: AtomMapper, scorer, n_centers: int =2,
                 n_processes: int = 1,
                 _initial_edge_lister: NetworkGenerator = None):
        """
        The Twin Star Ligand Network Planner , set's n ligands ligand into the center of a graph and connects all other ligands to each center.

        Parameters
        ----------
        mapper : AtomMapper
            the atom mapper is required, to define the connection between two ligands.
        scorer : AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        n_centers: int, optional
            the number of centers in the network. (default: 2)
        n_processes: int, optional
            number of processes that can be used for the network generation. (default: 1)
        _initial_edge_lister: LigandNetworkPlanner, optional
            this LigandNetworkPlanner is used to give the initial set of edges. For standard usage, the Maximal NetworPlanner is used.
            However in large scale approaches, it might be interesting to use the heuristicMaximalNetworkPlanner.. (default: MaximalNetworkPlanner)
        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(mapper=mapper,
                                                           scorer=scorer,
                                                           n_processes=n_processes)

        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=RadialNetworkAlgorithm(n_centers=n_centers),
                         n_processes=n_processes,
                         _initial_edge_lister=_initial_edge_lister)

        self.n_centers = n_centers


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
            components=components)
        mappings = initial_network.edges

        # Translate Mappings to graphable:
        edge_map = {(components.index(m.componentA),
                     components.index(m.componentB)): m for m in mappings}
        edges = list(sorted(edge_map.keys()))
        weights = [edge_map[k].annotations['score'] for k in edges]

        rg = self.network_generator.generate_network(edges=edges,
                                                     weights=weights)
        selected_mappings = [edge_map[k] for k in rg.edges]


        return LigandNetwork(edges=selected_mappings, nodes=components)
