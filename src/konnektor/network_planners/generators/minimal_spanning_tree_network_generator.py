# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from collections.abc import Iterable

from gufe import AtomMapper, Component, LigandNetwork

from konnektor.network_planners._networkx_implementations import MstNetworkAlgorithm

from ._abstract_network_generator import NetworkGenerator
from .maximal_network_generator import MaximalNetworkGenerator


class MinimalSpanningTreeNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
        scorer,
        n_processes: int = 1,
        progress: bool = False,
        _initial_edge_lister: NetworkGenerator = None,
    ):
        """
        The ``MinimalSpanningTreeNetworkGenerator``, builds an minimal spanning tree (MST) network for a given set of ``Component``/s.\
        The ``Transformation`` s of the network are represented by an ``AtomMapping`` s, which are scored by a ``AtomMappingScorer``.

        For the MST algorithm, the Kruskal Algorithm is used.

        The MST algorithm gives the optimal graph score possible and the minimal required set of ``Transformations``.
        This makes the  MST Network very efficient. However, the MST is not very robust, in case of one failing
        ``Transformation``, the network is immediately disconnected.
        The disconnectivity will translate to a loss of ``Component``/s in the final FE Network.

        Parameters
        ----------
        mapper :  Union[AtomMapper, list[AtomMapper]]
            ``AtomMapper`` or list of ``AtomMapper``/s to use to define the relationship between two ligands.
        scorer : AtomMappingScorer
            The scoring function to use for evaluating an atom mapping. Should give a score in [0,1].
        n_processes: int, optional
            Number of processes to be used for parallelization. (default: 1)
        progress: bool, optional
            If True, a progress bar will be displayed. (default: False)
        _initial_edge_lister: NetworkPlanner, optional
            ``NetworkPlanner`` to be used to generate the initial set of edges. For standard usage, the Maximal NetworkPlanner is often appropriate.
            For very large networks, the ``HeuristicMaximalNetworkPlanner`` might be a useful alternative.
            (default: MaximalNetworkPlanner)
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

    def generate_ligand_network(self, components: Iterable[Component]) -> LigandNetwork:
        """
        Generate a MST network from the given ``Component``/s.

        Parameters
        ----------
        components: Iterable[Component]
            ``Components`` to be used as nodes in the ``LigandNetwork``.

        Returns
        -------
        LigandNetwork
            ``LigandNetwork`` generated following the MST rules.

        """

        initial_network = self._initial_edge_lister.generate_ligand_network(components=components)
        mappings = initial_network.graph.edges(data=True)

        # Translate Mappings to graphable:
        edge_map = {
            (components.index(componentA), components.index(componentB)): d["object"]
            for componentA, componentB, d in mappings
        }
        edges = list(edge_map.keys())
        weights = [edge_map[k].annotations["score"] for k in edges]

        mg = self.network_generator.generate_network(edges, weights)

        # TODO: collect all the mappings, use j->i mapping if i->j not found? - double check this
        selected_mappings = [
            edge_map[k] if (k in edge_map) else edge_map[tuple(list(k)[::-1])] for k in mg.edges
        ]

        # intentionally make the ligand_network based *only* on the edges,
        # so we can catch any missing nodes in the next step
        mst_ligand_network = LigandNetwork(edges=selected_mappings)

        # check for a disconnected network
        missing_nodes = set(initial_network.nodes) - set(mst_ligand_network.nodes)
        if missing_nodes:
            raise RuntimeError(
                "ERROR: Unable to create edges for the following nodes: " + str(list(missing_nodes))
            )

        return mst_ligand_network
