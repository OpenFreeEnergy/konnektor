import itertools

from typing import Iterable, Callable, Union, Tuple
import functools
from tqdm.auto import tqdm

from gufe import SmallMoleculeComponent, AtomMapper

from openfe.setup.ligand_network import LigandNetwork    # only temproary
from ._abstract_ligand_network_planner import easyLigandNetworkPlanner


class ExplicitNetwork(easyLigandNetworkPlanner):
    def __init__(self, mapper, scorer=None):
        super().__init__(mapper=mapper, scorer=scorer,
                       network_generator=None)

    def generate_ligand_network(self,
                         edges: Iterable[Tuple[SmallMoleculeComponent,SmallMoleculeComponent]],
                         progress: Union[bool, Callable[[Iterable], Iterable]] = True,
                         ):
        """Create a network with pre-defined edges.

        """
        nodes = list(set([n for e in edges for n in e]))

        if progress is True:
            # default is a tqdm progress bar
            progress = functools.partial(tqdm, total=len(edges), delay=1.5)
        elif progress is False:
            progress = lambda x: x

        mappings = []
        for compoundA, compoundB in progress(edges):
            mappings = self.mapper.suggest_mappings(compoundA, compoundB)

        if self.scorer:
            mappings = [mapping.with_annotations({'score': self.scorer(mapping)})
                        for mapping in mappings]

        network = LigandNetwork(edges=mappings, nodes=nodes)
        return network