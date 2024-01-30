import itertools
import functools
import multiprocessing as mult

from tqdm.auto import tqdm

from typing import Iterable, Callable, Union, Optional

from gufe import SmallMoleculeComponent, AtomMapper

from konnektor.utils  import LigandNetwork    # only temproary
from ._abstract_ligand_network_planner import LigandNetworkPlanner


def thread_mapping(args):
    '''
    Helper function working as thread for parallel execution.
    Parameters
    ----------
    compound_pair

    Returns
    -------

    '''
    jobID, compound_pairs, mapper, scorer = args
    mapping_generator = [next(mapper.suggest_mappings(
        compound_pair[0], compound_pair[1])) for
        compound_pair in compound_pairs]

    if scorer:
        mappings = [mapping.with_annotations(
            {'score': scorer(mapping)})
            for mapping in mapping_generator]
    else:
        mappings = list(mapping_generator)

    return mappings

class MaximalNetworkPlanner(LigandNetworkPlanner):
    def __init__(self, mapper, scorer, progress=False, nprocesses=1):
        super().__init__(mapper=mapper, scorer=scorer,
                       network_generator=None, _initial_edge_lister=self)
        self.progress = progress
        self.nprocesses = nprocesses

    def generate_ligand_network(self,  nodes: Iterable[SmallMoleculeComponent],

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
        total = len(nodes) * (len(nodes) - 1) // 2
        n_batches = 1
        if self.progress is True and self.nprocesses<2:
            # default is a tqdm progress bar
            progress = functools.partial(tqdm, total=total, delay=1.5,
                                         desc="Mapping")
        elif self.progress is True and self.nprocesses>1:
            n_batches = 10*self.nprocesses
            progress = functools.partial(tqdm, total=n_batches, delay=1.5,
                                         desc="Mapping")
        else:
            progress = lambda x: x

        # Parallel or not Parallel:
        if(self.nprocesses > 1):
            # size of each batch +fetch division rest
            batch_num = (total//n_batches)+1

            # Prepare parallel execution.
            combinations =  list(itertools.combinations(nodes,2))
            batches = (combinations[i:i + n] for i in range(0, len(combinations), n_batches))

            jobs = [(job_id, combination, self.mapper, self.scorer) for job_id,
                    combination in enumerate(batches)]

            #Execute parallelism
            mappings = []
            with mult.Pool(self.nprocesses) as p:
                for sub_result in progress(p.imap(thread_mapping, jobs)):
                    mappings.extend(sub_result)
        else:
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
