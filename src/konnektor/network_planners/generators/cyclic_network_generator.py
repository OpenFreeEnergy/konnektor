# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from typing import Union, List, Iterable

from gufe import Component, LigandNetwork, AtomMapper

from konnektor.network_planners._networkx_implementations import \
    CyclicNetworkAlgorithm as nx_CNG
from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator


# Todo: check this algorithm again

class CyclicNetworkGenerator(NetworkGenerator):

    def __init__(self, mapper: AtomMapper, scorer,
                 node_present_in_cycles: int = 2,
                 cycle_sizes: Union[int, List[int]] = 3,
                 n_processes: int = 1,
                 _initial_edge_lister: NetworkGenerator = None):
        """the cyclic ligand planner tries to build up a network in which each
        node is contained in n cycles of a given size or size range.
        In order to do so, and to be time efficient, the class uses greedy
        algorithms to solve the problem.

        Parameters
        ----------
        mapper : AtomMapper
            the AtomMappers to use to propose mappings.  At least 1 required,
            but many can be given, in which case all will be tried to find the
            lowest score edges
        scorer : AtomMappingScorer
            any callable which takes a AtomMapping and returns a float
        node_present_in_cycles: int
            the number of cycles, the node should be present in.
        cycle_sizes: Union[int, List[int]]
            the cycle size in the graph, that is used for designing the graph.
            When providing a list[int], a range of sizes is allowed.
        n_processes: int, optional
            number of processes that can be used for the network generation.
            (default: 1)
        _initial_edge_lister: LigandNetworkPlanner, optional
            this LigandNetworkPlanner is used to give the initial set of edges.
            For standard usage, the Maximal NetworPlanner is used.
            However in large scale approaches, it might be interesting to use
            the heuristicMaximalNetworkPlanner. (default: MaximalNetworkPlanner)
        """

        network_generator = nx_CNG(
            node_cycle_connectivity=node_present_in_cycles,
            sub_cycle_size_range=cycle_sizes)
        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(mapper=mapper,
                                                           scorer=scorer,
                                                           n_processes=n_processes)

        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=network_generator,
                         n_processes=n_processes,
                         _initial_edge_lister=_initial_edge_lister)

    def generate_ligand_network(self, components: Iterable[
        Component]) -> LigandNetwork:
        """
           generate a cyclic network for the given compounds.


            Parameters
            ----------
            components: Iterable[Component]
                the components to be used for the LigandNetwork

            Returns
            -------
            LigandNetwork
                a cyclic network.
        """

        # Build Full Graph
        initial_networks = self._initial_edge_lister.generate_ligand_network(
            components=components)
        mappings = initial_networks.edges

        # Translate Mappings to graphable:
        # print("prepare network")
        edge_map = {
            (components.index(m.componentA), components.index(m.componentB)): m
            for m in mappings}
        edges = list(sorted(edge_map.keys()))
        weights = [edge_map[k].annotations['score'] for k in edges]

        # print("calculate Network")
        cg = self.network_generator.generate_network(edges=edges,
                                                     weights=weights)

        selected_mappings = [edge_map[k] for k in cg.edges]
        # print("Done")

        return LigandNetwork(edges=selected_mappings, nodes=components)
