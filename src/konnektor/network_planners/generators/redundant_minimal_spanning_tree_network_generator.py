# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from typing import Iterable

from gufe import Component, LigandNetwork, AtomMapper

from konnektor.network_planners._networkx_implementations import \
    MstNetworkAlgorithm
from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator


class RedundantMinimalSpanningTreeNetworkGenerator(NetworkGenerator):

    def __init__(self, mapper: AtomMapper, scorer, n_redundancy: int = 2,
                 n_processes: int = 1,
                 _initial_edge_lister: NetworkGenerator = None):
        """
        The `RedundantMinimalSpanningTreeNetworkGenerator` is an approach, that tries to increase the robustness over `Transformation` failures in an MST like network.
         
        The algorithm executes the MST algorithm `n_redundancy` times, always removes already selected `Transformations` in each iteration, and finally builds the overlay of all the newtorks. 
        This is constructing the Redundant MST Network. 
        
        In this way the number of edges is increased, but also the network is less vulnerable to `Transformation` failures.

        Parameters
        ----------
        mapper : AtomMapper
            the `AtomMapper`s to use to propose `AtomMapping`s.  At least 1 required,
            but many can be given, in which case all will be tried to find the
            lowest score edges
        scorer : AtomMappingScorer
            any callable which takes a `AtomMapping` and returns a float
        n_redundancy: int
            use MST n times to get a redundant set of `Transformations`.
        n_processes: int, optional
            number of processes that can be used for the network generation. (default: 1)
        _initial_edge_lister: NetworkGenerator, optional
            this `NetworkGenerator` is used to give the initial set of edges. For standard usage, the `MaximalNetworkGenerator` is used. (default: MaximalNetworkPlanner)


        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(mapper=mapper,
                                                           scorer=scorer,
                                                           n_processes=n_processes)

        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=MstNetworkAlgorithm(),
                         n_processes=n_processes,
                         _initial_edge_lister=_initial_edge_lister)

        self.n_redundancy = n_redundancy

    def generate_ligand_network(self, components: Iterable[
        Component]) -> LigandNetwork:
        """
        generate a redundant mst network for the given compounds.


        Parameters
        ----------
        components: Iterable[Component]
            the components to be used for the LigandNetwork

        Returns
        -------
        LigandNetwork
            a redundant MST network.
        """

        initial_network = self._initial_edge_lister.generate_ligand_network(
            components=components)
        mappings = initial_network.edges

        # Translate Mappings to graphable:
        edge_map = {
            (components.index(m.componentA), components.index(m.componentB)): m
            for m in mappings}
        edges = list(edge_map.keys())
        weights = [edge_map[k].annotations['score'] for k in edges]

        edge_weight = dict(zip(edges, weights))

        selected_edges = []
        for _ in range(self.n_redundancy):
            edges = list(edge_weight.keys())

            # filter for already selected edges
            filter_sEdges = lambda x: (x not in selected_edges and
                                       x[::-1] not in selected_edges)
            edges = list(filter(filter_sEdges, edges))
            # print("Edges", len(edges), edges)

            weights = [edge_weight[e] for e in edges]
            edge_weight = dict(list(zip(edges, weights)))
            # print("dEdges", len(edge_weight))
            ns = set([n for e in edges for n in e])
            # print("nodes", len(ns), ns)

            mg = self.network_generator.generate_network(edges, weights)

            if not mg.connected:
                nodes_index = {l: components.index(l) for l in components}
                missing_nodes = [l for l in components if
                                 (nodes_index[l] in mg.nodes)]
                raise RuntimeError("Unable to create edges for some nodes: "
                                   + str(list(missing_nodes)))
            # print("sel", len(mg.edges), mg.edges)

            selected_edges.extend(list(mg.edges))
            # print("collected", len(selected_edges), selected_edges)

        selected_mappings = [edge_map[k] if (k in edge_map) else edge_map[
            tuple(list(k)[::-1])] for k in selected_edges]

        return LigandNetwork(edges=selected_mappings, nodes=components)
