import logging
from typing import Iterable


from gufe import AtomMapper, AtomMappingScorer, LigandNetwork
from ._abstract_network_concatenator import NetworkConcatenator
from .._networkx_implementations import MstNetworkGenerator

log = logging.getLogger(__name__)


# Todo: check this algorithm again

class MstConcatenate(NetworkConcatenator):
    def __init__(self, mapper: AtomMapper, scorer: AtomMappingScorer, n_connecting_edges: int = 3,
                 node_present_in_cycles: int = 2, cycle_sizes: Union[int, List[int]] = 3, nprocesses: int = 1):
        """
        This concatenators is connnecting two Networks with a kruskal like approach up to the number of connecting edges.

        Parameters
        ----------
        mapper: AtomMapper
            the atom mapper is required, to define the connection between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        n_connecting_edges: int, optional
            number of connecting edges. (default: 3)
        nprocesses: int
            number of processes that can be used for the network generation. (default: 1)
        """
        super().__init__(mapper=mapper, scorer=scorer, network_generator=MstNetworkGenerator(), nprocesses=nprocesses)
        self.n_connecting_edges = n_connecting_edges

    def concatenate_networks(self, ligand_networks: Iterable[LigandNetwork]) -> LigandNetwork:
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
            returns a concatenated LigandNetwork object, containing all networks.

        """
        raise NotImplementedError()

        # Todo: implement.
