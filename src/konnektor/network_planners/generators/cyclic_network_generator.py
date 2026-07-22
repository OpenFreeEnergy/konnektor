# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from collections.abc import Callable, Iterable

from gufe import AtomMapper, AtomMapping, Component, LigandNetwork

from konnektor.network_planners._networkx_implementations import (
    CyclicNetworkAlgorithm as nx_CNG,
)

from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator

# Todo: check this algorithm again


class CyclicNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
        scorer: Callable[[AtomMapping], float] | None,
        node_present_in_cycles: int = 2,
        cycle_sizes: int | list[int] = 3,
        n_processes: int = 1,
        progress: bool = False,
        _initial_edge_lister: NetworkGenerator | None = None,
    ):
        """
        A `NetworkGenerator` that generates a network containing cycles.

        The greedy algorithm builds the network up from a node-wise perspective.
        For each node, the algorithm generates all cycles of size `cycle_sizes` and assigns a score to each cycle as the sum of all sub-scores.
        Next, it selects the `node_present_in_cycles` best score performing and node diversity increasing (see below) cycles per node.
        This set of selected edges are used to construct the graph.
        The node diversity criterion is an addition which biases to spread the cycles on the graph equally between all `Components`.

        The number of cycles around each `Component` can be defined by `component_present_in_cycles`
        and the allowed cycle size can be tweaked with `cycle_sizes`.

        This layout has well-distributed connectivity between all `Component` s, which increases the robustness very well,
        but still allows for a better graph score then the Twin Star Network, as the connectivity distribution is biased and not enforced.
        The large number of cycles might be very useful for statistical analysis, but does come with the added cost of additional edges.


        Parameters
        ----------
        mappers : AtomMapper | list[AtomMapper]
            AtomMapper(s) to use to propose mappings. If more than one AtomMapper is provided, all will be tried to find the
            lowest score for each edges.
        scorer :  Callable[[AtomMapping], float] | None
            Scoring function that takes an AtomMapping and returns a score in [0,1].
        node_present_in_cycles: int
            Number of cycles each node should be present in.
        cycle_sizes: int | list[int]
            The cycle size to be used for designing the graph.
            When providing a list[int], a range of sizes is allowed (e.g. `[3,4]`). (default: 3)
        n_processes: int, optional
            Number of processes that can be used for the network generation, default 1.
        progress: bool, optional
            If True, displays a progress bar, default False.
        _initial_edge_lister: NetworkGenerator | None
            The NetworkGenerator to use  if the NetworkGenerator requires an initial set of edges.
            If None, a ``MaximalNetworkGenerator`` will be used with the provided mappers and scorer.
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
        """Generate a cyclic network for the given compounds.

        Parameters
        ----------
        components : Iterable[Component]
            ``Component``/s to include as nodes in the ``LigandNetwork``.

        Returns
        -------
        LigandNetwork
            The generated cyclic ``LigandNetwork``.

        """

        # build initial graph
        initial_networks = self._initial_edge_lister.generate_ligand_network(components=components)
        mappings = initial_networks.edges

        # Translate Mappings to graphable:
        edge_map = {
            (components.index(m.componentA), components.index(m.componentB)): m for m in mappings
        }
        edges = list(sorted(edge_map.keys()))
        weights = [edge_map[k].annotations["score"] for k in edges]

        cg = self.network_generator.generate_network(edges=edges, weights=weights)

        selected_mappings = [edge_map[k] for k in cg.edges]

        return LigandNetwork(edges=selected_mappings, nodes=components)
