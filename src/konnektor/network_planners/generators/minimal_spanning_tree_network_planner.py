from typing import Iterable

from gufe import LigandNetwork, Component, AtomMapper, AtomMappingScorer
from konnektor.network_planners._networkx_implementations import MstNetworkGenerator

from ._abstract_network_generator import NetworkGenerator
from .maximal_network_planner import MaximalNetworkGenerator


class MinimalSpanningTreeNetworkGenerator(NetworkGenerator):

    def __init__(self, mapper: AtomMapper, scorer: AtomMappingScorer,
                 nprocesses: int = 1, _initial_edge_lister: NetworkGenerator = None):
        """
        The minimal spanning tree ligand network planner, builds an MST for a given set of ligands. The edges of the the graph,
        are represented by an AtomMapping, which is scored by the AtomMappingScorer. The MST algorithm gives in theory the most efficient graph possible.
        However, the MST is not very robust, in case of one failing edge.

        Parameters
        ----------
        mapper : AtomMapper
            the atom mapper is required, to define the connection between two ligands.
        scorer : AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        nprocesses: int, optional
            number of processes that can be used for the network generation. (default: 1)
        _initial_edge_lister: LigandNetworkPlanner, optional
            this LigandNetworkPlanner is used to give the initial set of edges. For standard usage, the Maximal NetworPlanner is used.
            However in large scale approaches, it might be interesting to use the heuristicMaximalNetworkPlanner.. (default: MaximalNetworkPlanner)
        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(mapper=mapper, scorer=scorer, nprocesses=nprocesses)

        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=MstNetworkGenerator(),
                         nprocesses=nprocesses,
                         _initial_edge_lister=_initial_edge_lister)

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
        """
        Generate a MST network from a list of components.

        Parameters
        ----------
        components: Iterable[Component]
        the components to be used for the LigandNetwork

        Returns
        -------
        LigandNetwork
            a ligand network following the MST rules.

        """

        initial_network = self._initial_edge_lister.generate_ligand_network(
            components=components)
        mappings = initial_network.edges

        # Translate Mappings to graphable:
        edge_map = {(components.index(m.componentA), components.index(m.componentB)): m for m in mappings}
        edges = list(edge_map.keys())
        weights = [edge_map[k].annotations['score'] for k in edges]

        mg = self.network_generator.generate_network(edges, weights)

        if not mg.connected:
            nodes_index = {l: components.index(l) for l in components}
            missing_nodes = [l for l in components if (nodes_index[l] in mg.nodes)]
            raise RuntimeError("Unable to create edges for some nodes: "
                               + str(list(missing_nodes)))

        selected_mappings = [edge_map[k] if (k in edge_map) else edge_map[tuple(list(k)[::-1])] for k in mg.edges]

        return LigandNetwork(edges=selected_mappings, nodes=components)
