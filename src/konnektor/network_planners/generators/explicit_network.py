import functools
from typing import Iterable, Tuple

from gufe import Component, LigandNetwork, AtomMapper, AtomMappingScorer
from tqdm.auto import tqdm

from ._abstract_ligand_network_generator import LigandNetworkGenerator
from ._parallel_mapping_pattern import _parallel_map_scoring


class ExplicitNetwork(LigandNetworkGenerator):
    def __init__(self, mapper: AtomMapper, scorer: AtomMappingScorer, progress: bool = False, ):
        """

        Parameters
        ----------
        mapper: AtomMapper
            the atom mapper is required, to define the connection between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        """
        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=None)
        self.progress = progress

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

        if self.nprocesses > 1:  # go parallel
            mappings = _parallel_map_scoring(
                possible_edges=edges,
                scorer=self.scorer,
                mapper=self.mapper,
                n_processes=self.nprocesses,
                show_progress=self.progress)
        else:  # sequential implementation
            if self.progress:
                # default is a tqdm progress bar
                progress = functools.partial(tqdm, total=len(edges), delay=1.5)
            else:
                progress = lambda x: x

            mappings = []
            for compoundA, compoundB in progress(edges):
                mappings = self.mapper.suggest_mappings(compoundA, compoundB)

            if self.scorer:
                mappings = [mapping.with_annotations({'score': self.scorer(mapping)})
                            for mapping in mappings]

        network = LigandNetwork(edges=mappings, nodes=nodes)
        return network
