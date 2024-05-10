from typing import Iterable, Tuple

from gufe import Component, LigandNetwork, AtomMapper, AtomMappingScorer

from ._abstract_network_generator import NetworkGenerator
from ._parallel_mapping_pattern import _parallel_map_scoring


class ExplicitNetworkGenerator(NetworkGenerator):
    def __init__(self, mapper: AtomMapper, scorer: AtomMappingScorer,
                 n_processes: int = 1, show_progress: bool = False, ):
        """

        Parameters
        ----------
        mapper: AtomMapper
            the atom mapper is required, to define the connection between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        n_processes: int
            number of processes used to build the ligand network
        show_progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        """
        super().__init__(mapper=mapper, scorer=scorer, nprocesses=n_processes,
                         network_generator=None)
        self.progress = show_progress

    def generate_ligand_network(self,
                                edges: Iterable[Tuple[Component, Component]],
                                ) -> LigandNetwork:
        """
        Create a network with pre-defined edges.

        This class is can be used as initial_edge_lister

        Parameters
        ----------
        edges: Iterable[Tuple[Component,Component]]
            planned edges, that will be connected with mappings and scores.
            Each Tuple in this case represent one edge.

        Returns
        -------
        LigandNetwork
            the provided network.
        """
        nodes = list(set([n for e in edges for n in e]))

        mappings = _parallel_map_scoring(
            possible_edges=edges,
            scorer=self.scorer,
            mapper=self.mapper,
            n_processes=self.nprocesses,
            show_progress=self.progress)

        network = LigandNetwork(edges=mappings, nodes=nodes)
        return network
