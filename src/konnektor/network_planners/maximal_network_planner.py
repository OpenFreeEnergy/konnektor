import itertools

from typing import Iterable, Callable, Union, Optional
import functools
from tqdm.auto import tqdm

from gufe import SmallMoleculeComponent, AtomMapper

from konnektor.utils  import LigandNetwork    # only temproary
from ._abstract_ligand_network_planner import easyLigandNetworkPlanner


class MaximalNetworkPlanner(easyLigandNetworkPlanner):
    def __init__(self, mapper, scorer):
        super().__init__(mapper=mapper, scorer=scorer,
                       network_generator=None)

    def generate_ligand_network(self,
                         nodes: Iterable[SmallMoleculeComponent],
                         progress: Union[bool, Callable[[Iterable], Iterable]] = True,
                         # allow_disconnected=True
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
        nodes = list(nodes)

        if progress is True:
            # default is a tqdm progress bar
            total = len(nodes) * (len(nodes) - 1) // 2
            progress = functools.partial(tqdm, total=total, delay=1.5)
        elif progress is False:
            progress = lambda x: x
        # otherwise, it should be a user-defined callable

        ligands, mappings = self._input_generate_all_possible_mappings(ligands=nodes)

        mapping_generator = itertools.chain.from_iterable(
            self.mapper.suggest_mappings(molA, molB)
            for molA, molB in progress(itertools.combinations(nodes, 2))
        )
        if self.scorer:
            mappings = [mapping.with_annotations({'score': self.scorer(mapping)})
                        for mapping in mapping_generator]
        else:
            mappings = list(mapping_generator)

        network = LigandNetwork(mappings, nodes=nodes)
        return network