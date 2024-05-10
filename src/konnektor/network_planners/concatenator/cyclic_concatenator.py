import logging
from typing import Iterable

from gufe import AtomMapper, AtomMappingScorer, LigandNetwork

from ._abstract_network_concatenator import NetworkConcatenator
from .max_concatenator import MaxConcatenate
from .._networkx_implementations import MstNetworkGenerator

log = logging.getLogger(__name__)


# Todo: check this algorithm again

class CyclicConcatenate(NetworkConcatenator):
    def __init__(self, mapper: AtomMapper, scorer: AtomMappingScorer,
                 n_connecting_cycles: int = 2,
                 cycle_sizes: Union[int, List[int]] = 3, nprocesses: int = 1,
                 _initial_edge_lister: NetworkConcatenator = None):
        """
        This concatenator is connnecting two Networks with a kruskal like
        approach up to the number of connecting edges.

        Parameters
        ----------
        mapper: AtomMapper
            the atom mapper is required, to define the connection
            between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a
            score between [0,1].
        n_connecting_cycles: int, optional
            build at least n cycles between th networks. (default: 2)
        cycle_sizes: Union[int, List[int]], optional
            build cycles of given size. (default:3)
        nprocesses: int
            number of processes that can be used for the network generation.
            (default: 1)
        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaxConcatenate(mapper=mapper, scorer=scorer,
                                                  nprocesses=nprocesses)

        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=MstNetworkGenerator(),
                         nprocesses=nprocesses)
        self.n_connecting_edges = n_connecting_cycles
        self.cycle_sizes = cycle_sizes

    def concatenate_networks(self, ligand_networks: Iterable[
        LigandNetwork]) -> LigandNetwork:
        """

        Parameters
        ----------
        ligand_networks: Iterable[LigandNetwork]
            an iterable of ligand networks, that shall be connected.
        n_connecting_edges: int
            number of edges, to connect the networks

        Returns
        -------
        LigandNetwork
            returns a concatenated LigandNetwork object,
            containing all networks.

        """
        raise NotImplementedError()
        # Todo: implement.

        selected_edges = []
        selected_nodes = []
        for ligandNetworkA, ligandNetworkB in itertools.combinations(
                ligand_networks, 2):
            # Generate fully connected Bipartite Graph
            ligands = list(ligandNetworkA.nodes) + list(ligandNetworkB.nodes)
            fully_connected_graph = self._initial_edge_lister(
                [ligandNetworkA, ligandNetworkB])
            bipartite_graph_mappings = list(fully_connected_graph.edges)

            # TODO Cycle Selection

            selected_edges.extend(selected_mappings)

        # Constructed final Edges:
        # Add all old network edges:
        for network in ligand_networks:
            selected_edges.extend(network.edges)
            selected_nodes.extend(network.nodes)

        concat_LigandNetwork = LigandNetwork(edges=selected_edges,
                                             nodes=set(selected_nodes))

        log.info("Total Concatenated Edges:  " + str(len(selected_edges)))

        return concat_LigandNetwork
