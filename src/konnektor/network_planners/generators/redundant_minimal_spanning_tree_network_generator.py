# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import warnings
from collections.abc import Iterable

import networkx as nx
from gufe import AtomMapper, Component, LigandNetwork

from konnektor.network_planners._networkx_implementations import MstNetworkAlgorithm

from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator


class RedundantMinimalSpanningTreeNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
        scorer,
        n_redundancy: int = 2,
        n_processes: int = 1,
        progress: bool = False,
        _initial_edge_lister: NetworkGenerator = None,
    ):
        """
        The ``RedundantMinimalSpanningTreeNetworkGenerator`` is an approach that tries to increase
        robustness to ``Transformation`` failures in an MST-like network.

        This algorithm executes the MST algorithm ``n_redundancy`` times, always removing
        already-selected ``Transformations`` in each iteration, and finally builds the overlay of all the networks.
        This is constructs the Redundant MST Network.

        In this way, the number of edges is increased, but the network is also less vulnerable to ``Transformation`` failures.

        Parameters
        ----------
        mappers :  Union[AtomMapper, list[AtomMapper]]
            The ``AtomMapper``/s to use to propose ``AtomMapping``/s.  At least 1 required,
            but many can be given, in which case all will be tried to find the
            lowest score edges
        scorer : AtomMappingScorer
            any callable which takes a ``AtomMapping`` and returns a float
        n_redundancy: int
            use MST n times to get a redundant set of ``Transformations``.
        n_processes: int, optional
            number of processes that can be used for the network generation. (default: 1)
        progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        _initial_edge_lister: NetworkGenerator, optional
            this ``NetworkGenerator`` is used to give the initial set of edges. For standard usage, the ``MaximalNetworkGenerator`` is used. (default: MaximalNetworkPlanner)


        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(
                mappers=mappers, scorer=scorer, n_processes=n_processes
            )

        super().__init__(
            mappers=mappers,
            scorer=scorer,
            network_generator=MstNetworkAlgorithm(),
            n_processes=n_processes,
            progress=progress,
            _initial_edge_lister=_initial_edge_lister,
        )

        self.n_redundancy = n_redundancy

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
        """
        Generate a redundant MST network from the given ``Component``/s.


        Parameters
        ----------
        components: Iterable[Component]
            ``Components`` to be used as nodes in the ``LigandNetwork``.

        Returns
        -------
        LigandNetwork
            a redundant MST network.
        """

        initial_network = self._initial_edge_lister.generate_ligand_network(components=components)
        mappings = initial_network.edges

        # Translate Mappings to graphable:
        edge_map = {
            (components.index(m.componentA), components.index(m.componentB)): m for m in mappings
        }
        edges = list(edge_map.keys())
        weights = [edge_map[k].annotations["score"] for k in edges]

        edge_weight = dict(zip(edges, weights))

        selected_edges = []
        for n in range(self.n_redundancy):
            edges = list(edge_weight.keys())

            # filter for already selected edges
            filter_sEdges = lambda x: (x not in selected_edges and x[::-1] not in selected_edges)
            edges = list(filter(filter_sEdges, edges))

            weights = [edge_weight[e] for e in edges]
            edge_weight = dict(list(zip(edges, weights)))
            # TODO: if weights and edge_weight are [], this breaks!
            # try/except for now to reproduce openfe behavior, but handle this better
            try:
                mg = self.network_generator.generate_network(edges, weights)
            except nx.NetworkXPointlessConcept:
                warnings.warn(
                    f"Cannot create any minimal spanning network for redundancy iteration {n + 1}"
                )
            # TODO: im not sure this will ever catch isolated nodes, since it's constructed explicitly from edges
            if not mg.connected:
                nodes_index = {c: components.index(c) for c in components}
                missing_nodes = [c for c in components if (nodes_index[c] in mg.nodes)]
                raise RuntimeError(
                    f"Unable to create edges for the following node during redundancy iteration {n + 1}: {missing_nodes}"
                )

            selected_edges.extend(list(mg.edges))
        selected_mappings = [
            edge_map[k] if (k in edge_map) else edge_map[tuple(list(k)[::-1])]
            for k in selected_edges
        ]

        # intentionally make the ligand_network based *only* on the edges,
        # so we can catch any missing nodes in the next step
        ligand_network = LigandNetwork(edges=selected_mappings)

        # check for a disconnected network
        missing_nodes = set(initial_network.nodes) - set(ligand_network.nodes)
        if missing_nodes:
            raise RuntimeError(
                "ERROR: Unable to create edges for the following nodes: " + str(list(missing_nodes))
            )

        return ligand_network
