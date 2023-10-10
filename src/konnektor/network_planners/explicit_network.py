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
        """Create a network with all possible proposed mappings.

        This will attempt to create (and optionally score) all possible mappings
        (up to $N(N-1)/2$ for each mapper given). There may be fewer actual
        mappings that this because, when a mapper cannot return a mapping for a
        given pair, there is simply no suggested mapping for that pair.
        This network is typically used as the starting point for other network
        generators (which then optimize based on the scores) or to debug atom
        mappers (to see which mappings the mapper fails to generate).


        Parameters
        ----------
        nodes : Iterable[SmallMoleculeComponent]
          the ligands to include in the LigandNetwork
        mappers : Iterable[LigandAtomMapper]
          the AtomMappers to use to propose mappings.  At least 1 required,
          but many can be given, in which case all will be tried to find the
          lowest score edges
        scorer : Scoring function
          any callable which takes a LigandAtomMapping and returns a float
        progress : Union[bool, Callable[Iterable], Iterable]
          progress bar: if False, no progress bar will be shown. If True, use a
          tqdm progress bar that only appears after 1.5 seconds. You can also
          provide a custom progress bar wrapper as a callable.
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