# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from typing import Union, List, Iterable

from gufe import Component, LigandNetwork, AtomMapper

from konnektor.network_planners._networkx_implementations import (
    CyclicNetworkAlgorithm as nx_CNG,
)
from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator


# Todo: check this algorithm again


class CyclicNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: Union[AtomMapper, list[AtomMapper]],
        scorer,
        node_present_in_cycles: int = 2,
        cycle_sizes: Union[int, List[int]] = 3,
        n_processes: int = 1,
        progress: bool = False,
        _initial_edge_lister: NetworkGenerator = None,
    ):
        """
        A `NetworkGenerator` that generates a network based on many network cycles.
        This is of interest for analyzing the uncertainty of FE Estimates along thermodynamic
        cycles and possibly for correcting the estimates with cycle closure analysis.

        The greedy algorithm builds the network up from a nodewise perspective.
        For each node, the algorithm generates all cycles of size `cycle_sizes` and assigns a score to each cycle as the sum of all sub-scores.
        Next, it selects the `node_present_in_cycles` best score performing and node diversity increasing (see below) cycles per node.
        The set of selected `Transformations` constructs the graph.
        The node diversity criterion is an addition which biases to spread the cycles on the graph equally between all `Components`.

        The number of cycles around each `Component` can be defined by `component_present_in_cycles`
        and the allowed cycle size can be tweaked with `cycle_sizes`.

        This layout has well-distributed connectivity between all `Component` s, which increases the robustness very well,
        but still allows for a better graph score then the Twin Star Network, as the connectivity distribution is biased not enforced.
        The large number of cycles might be very useful for statistical analysis.  Nevertheless, the network has an increased amount of `Transformation`s.


        Parameters
        ----------
        mappers : Union[AtomMapper, list[AtomMapper]]
            AtomMapper(s) to use to propose mappings.  At least 1 required,
            but many can be given, in which case all will be tried to find the
            lowest score edges.
        scorer : AtomMappingScorer
            Any callable which takes a AtomMapping and returns a float.
        node_present_in_cycles: int
            The number of cycles the node should be present in.
        cycle_sizes: Union[int, List[int]]
            The cycle size to be used for designing the graph.
            When providing a list[int], a range of sizes is allowed (e.g. `[3,4]`). (default: 3)
        n_processes: int, optional
            Number of processes that can be used for the network generation.
            (default: 1)
        progress: bool, optional
            If `True`, displays a progress bar. (default: False)
        _initial_edge_lister: NetworkPlanner, optional
            The NetworkPlanner used to give the initial set of edges.
            For standard usage, the MaximalNetworkPlanner is used. (default: MaximalNetworkPlanner)
        """

        network_generator = nx_CNG(
            node_cycle_connectivity=node_present_in_cycles,
            sub_cycle_size_range=cycle_sizes,
        )

        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(
                mappers=mappers, scorer=scorer, n_processes=n_processes
            )

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=network_generator,
            n_processes=n_processes,
            progress=progress,
            _initial_edge_lister=_initial_edge_lister,
        )

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
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
            components=components
        )
        mappings = initial_networks.edges

        # Translate Mappings to graphable:
        # print("prepare network")
        edge_map = {
            (components.index(m.componentA), components.index(m.componentB)): m
            for m in mappings
        }
        edges = list(sorted(edge_map.keys()))
        weights = [edge_map[k].annotations["score"] for k in edges]

        # print("calculate Network")
        cg = self.network_generator.generate_network(edges=edges, weights=weights)

        # cg = self.network_generator.generate_network_double_greedy(edges=edges,
        #                                                           weights=weights)

        selected_mappings = [edge_map[k] for k in cg.edges]
        # print("Done")

        return LigandNetwork(edges=selected_mappings, nodes=components)
