# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from typing import Iterable, Union

from gufe import Component, LigandNetwork, AtomMapper

from konnektor.network_planners._networkx_implementations import (
    NNodeEdgesNetworkAlgorithm,
)
from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator


class NNodeEdgesNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: Union[AtomMapper, list[AtomMapper]],
        scorer,
        target_component_connectivity: int = 3,
        n_processes: int = 1,
        progress: bool = False,
        _initial_edge_lister: NetworkGenerator = None,
    ):
        """
        The N-Node Edges Network tries to add more redundancy to the MST Network and tries to improve the robustness.

        The algorithm first build a MST Network. Then it will add best score performing `Transformations`
        to guarantee a 'Component' connectivity of `target_node_connectivity`.


        Parameters
        ----------
        mapper : AtomMapper
            the `AtomMapper` to use to propose `AtomMapping` s.
        scorer : AtomMappingScorer
            any callable which takes a `AtomMapping` and returns a float
        target_node_connectivity: int
            the number of connecting `Transformations` per `Component`.
        n_processes: int, optional
            number of processes that can be used for the network generation. (default: 1)
        progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        _initial_edge_lister: LigandNetworNetworkGeneratorkPlanner, optional
            this `NetworkGenerator` is used to give the initial set of `Transformation` s.
            For standard usage, the MaximalNetworGenerator is used, which will provide all possible `Transformation` s. (default: MaximalNetworkPlanner)

        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(
                mappers=mappers, scorer=scorer, n_processes=n_processes
            )

        network_generator = NNodeEdgesNetworkAlgorithm(
            target_node_connectivity=target_component_connectivity
        )
        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=network_generator,
            n_processes=n_processes,
            progress=progress,
            _initial_edge_lister=_initial_edge_lister,
        )

    @property
    def target_node_connectivity(self) -> int:
        """
        this property defines how many edges should be created per node.

        Returns
        -------
        int
            the targeted number of edges per node
        """
        return self.network_generator.target_node_connectivity

    @target_node_connectivity.setter
    def target_node_connectivity(self, target_node_connectivity: int):
        self.network_generator.target_node_connectivity = target_node_connectivity

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
        """Plan a Network which connects all ligands with at least n edges

        Parameters
        ----------
        components : Iterable[Component]
        the ligands to include in the Network

        Returns
        -------
        LigandNetwork
            the resulting ligand network.
        """
        # Build Full Graph
        initial_networks = self._initial_edge_lister.generate_ligand_network(
            components=components
        )
        mappings = initial_networks.edges

        # Translate Mappings to graphable:
        edge_map = {
            (components.index(m.componentA), components.index(m.componentB)): m
            for m in mappings
        }
        edges = list(sorted(edge_map.keys()))
        weights = [edge_map[k].annotations["score"] for k in edges]

        sg = self.network_generator.generate_network(edges=edges, weights=weights)

        selected_mappings = [
            edge_map[k] if (k in edge_map) else edge_map[tuple(list(k)[::-1])]
            for k in sg.edges
        ]
        return LigandNetwork(edges=selected_mappings, nodes=components)
