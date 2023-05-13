import itertools

from typing import Iterable, Callable, Union, Optional

from gufe import SmallMoleculeComponent, AtomMapper

from openfe.setup.ligand_network import LigandNetwork    # only temproary
from ._abstract_network_planner import _AbstractNetworkPlanner, Network


class maximalNetworkPlanner(_AbstractNetworkPlanner):
    def generate_network(
        ligands: Iterable[SmallMoleculeComponent],
        mappers: Iterable[AtomMapper],
        scorer: Optional[Callable[[AtomMapper], float]] = None,
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
        ligands : Iterable[SmallMoleculeComponent]
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
        nodes = list(ligands)

        if progress is True:
            # default is a tqdm progress bar
            total = len(nodes) * (len(nodes) - 1) // 2
            progress = functools.partial(tqdm, total=total, delay=1.5)
        elif progress is False:
            progress = lambda x: x
        # otherwise, it should be a user-defined callable

        mapping_generator = itertools.chain.from_iterable(
            mapper.suggest_mappings(molA, molB)
            for molA, molB in progress(itertools.combinations(nodes, 2))
            for mapper in mappers
        )
        if scorer:
            mappings = [mapping.with_annotations({'score': scorer(mapping)})
                        for mapping in mapping_generator]
        else:
            mappings = list(mapping_generator)

        network = LigandNetwork(mappings, nodes=nodes)
        return network