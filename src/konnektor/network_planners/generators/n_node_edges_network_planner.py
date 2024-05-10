from typing import Iterable

from gufe import Component, LigandNetwork, AtomMapper
from konnektor.network_planners._networkx_implementations import NNodeEdgesNetworkGenerator

from ._abstract_network_generator import NetworkGenerator
from .maximal_network_planner import MaximalNetworkGenerator


class NNodeEdgesNetworkGenerator(NetworkGenerator):

    def __init__(self, mapper: AtomMapper, scorer,
                 target_node_connectivity: int = 3,
                 nprocesses: int = 1,
                 _initial_edge_lister: NetworkGenerator = None):
        """
        the NNodeEdges Network planner, is building a graph, in which each node is connected by at least the target_node_connectivity.

        Parameters
        ----------
        mapper : AtomMapper
            the AtomMappers to use to propose mappings.  At least 1 required,
            but many can be given, in which case all will be tried to find the
            lowest score edges
        scorer : AtomMappingScorer
            any callable which takes a AtomMapping and returns a float
        target_node_connectivity: int
            the number of connecting edges per node.
        nprocesses: int, optional
            number of processes that can be used for the network generation. (default: 1)
        _initial_edge_lister: LigandNetworkPlanner, optional
            this LigandNetworkPlanner is used to give the initial set of edges. For standard usage, the Maximal NetworPlanner is used.
            However in large scale approaches, it might be interesting to use the heuristicMaximalNetworkPlanner.. (default: MaximalNetworkPlanner)

        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(mapper=mapper, scorer=scorer, nprocesses=nprocesses)

        network_generator = NNodeEdgesNetworkGenerator(target_node_connectivity=target_node_connectivity)
        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=network_generator,
                         nprocesses=nprocesses,
                         _initial_edge_lister=_initial_edge_lister)

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
            components=components)
        mappings = initial_networks.edges

        # Translate Mappings to graphable:
        edge_map = {(components.index(m.componentA), components.index(m.componentB)): m for m in mappings}
        edges = list(sorted(edge_map.keys()))
        weights = [edge_map[k].annotations['score'] for k in edges]

        sg = self.network_generator.generate_network(edges=edges, weights=weights)

        selected_mappings = [edge_map[k] if (k in edge_map) else edge_map[tuple(list(k)[::-1])] for k in sg.edges]
        return LigandNetwork(edges=selected_mappings, nodes=components)
