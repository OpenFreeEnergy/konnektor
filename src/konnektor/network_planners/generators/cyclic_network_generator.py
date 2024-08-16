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
        """
        The `CyclicNetworkGenerator` is generating a network based on many network cycles. T
        his is of interest for analyzing the uncertainty of FE Estimates along the thermodynamic cycles and possibly correct the estimates with cycle closure analysis.

        The greedy algorithm builds the network up from a nodewise perspective.
        For each node, the algorithm generates all cycles of size `cycle_sizes` and assigns a score to each cylce as the sum of all sub-scores. 
        Next it selects the `node_present_in_cyces` best score perfoming and node diversity increasing (see below) cycles per node.
        The set of selected `Transformations` constructs the graph.
        The node diversity criterion is an addition, which biases to spread the cycles on the graph eaqually between all `Components`
         
        The number of cylces, around each `Component` can be defined by `component_present_in_cycles` and allowed cylce size can be tweaked with `cycle_sizes`. For `cycle_sizes` either an integer for providing an expected cycle size (e.g. `3`) or a range of allowed cycle sizes (e.g. `[3,4]`).
        
        This layout has a well distributed connectivity between all `Component`s which increases the robustness very well, but still allows for a better graph score then the Twin Star Network, as the connectivity distribution is biased not enforced.
        The large number of cycles might be very useful for statical analysis.  Nevertheless, the network has an increased amount of `Transformation`s  
        

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
            When providing a list[int], a range of sizes is allowed (e.g. `[3,4]`). (default: 3)
        n_processes: int, optional
            number of processes that can be used for the network generation.
            (default: 1)
        _initial_edge_lister: LigandNetworkPlanner, optional
            this LigandNetworkPlanner is used to give the initial set of edges.
            For standard usage, the Maximal NetworPlanner is used. (default: MaximalNetworkPlanner)
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
        components : Iterable[Component]
          the ligands to include in the LigandNetwork

        Returns
        -------
        LigandNetwork
            a complex network.
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
        
        #cg = self.network_generator.generate_network_double_greedy(edges=edges, 
        #                                                           weights=weights)

        selected_mappings = [edge_map[k] for k in cg.edges]
        # print("Done")

        return LigandNetwork(edges=selected_mappings, nodes=components)
