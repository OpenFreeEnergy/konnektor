# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import functools
import multiprocessing as mult
from collections.abc import Callable

from gufe import AtomMapper, AtomMapping, SmallMoleculeComponent
from tqdm.auto import tqdm


def thread_mapping(args) -> list[AtomMapping]:
    """
    Helper function working as thread for parallel execution.

    Parameters
    ----------
    args:
        contains a list of: jobID, compound_pairs, mapper, scorer

    Returns
    -------
    list[AtomMapping]:
        return a list of scored atom mappings

    """

    jobID, compound_pairs, mappers, scorer = args  # noqa

    mappings = []
    for component_pair in compound_pairs:
        best_score = 0.0
        best_mapping = None
        molA = component_pair[0]
        molB = component_pair[1]

        for mapper in mappers:
            mapping_generator = mapper.suggest_mappings(molA, molB)

            if scorer:
                try:
                    tmp_mappings = [
                        mapping.with_annotations({"score": scorer(mapping)})
                        for mapping in mapping_generator
                    ]
                except:
                    continue
                if len(tmp_mappings) > 0:
                    tmp_best_mapping = min(tmp_mappings, key=lambda m: m.annotations["score"])

                    if tmp_best_mapping.annotations["score"] < best_score or best_mapping is None:
                        best_score = tmp_best_mapping.annotations["score"]
                        best_mapping = tmp_best_mapping

            else:
                try:
                    best_mapping = next(mapping_generator)
                except:
                    continue
        if best_mapping is not None:
            mappings.append(best_mapping)

    return mappings


def _parallel_map_scoring(
    possible_edges: list[tuple[SmallMoleculeComponent, SmallMoleculeComponent]],
    scorer: Callable[[AtomMapping], float],
    mappers: list[AtomMapper],
    n_processes: int,
    show_progress: bool = True,
) -> list[AtomMapping]:
    """
    This helper function parallelize mapping and scoring of a given list of
    molecule pairs.

    Parameters
    ----------
    possible_edges: tuple[SmallMoleculeComponent, SmallMoleculeComponent]
        two  molecules to be mapped.
    scorer: callable
        scoring the mappings
    mapper: AtomMapper
        atom mapper for the mappings
    n_processes: int
        number of processes for parallelization
    show_progress: bool
        show a tqdm progressbar.

    Returns
    -------
    list[AtomMapping]:
        return a list of scored atom mappings
    """
    if show_progress is True and n_processes > 1:
        n_batches = 10 * n_processes
        progress = functools.partial(tqdm, total=n_batches, delay=1.5, desc="Mapping")
    else:
        progress = lambda x: x

    possible_edges = list(possible_edges)
    n_batches = 10 * n_processes
    # total = len(possible_edges)

    # # size of each batch +fetch division rest
    # batch_num = (total // n_batches) + 1

    # Prepare parallel execution.
    # suboptimal implementation, but itertools.batch is python 3.12,
    batches = (possible_edges[i : i + n_batches] for i in range(0, len(possible_edges), n_batches))

    jobs = [(job_id, combination, mappers, scorer) for job_id, combination in enumerate(batches)]

    # Execute parallelism
    mappings = []
    with mult.Pool(n_processes) as p:
        for sub_result in progress(p.imap(thread_mapping, jobs)):
            mappings.extend(sub_result)

    return mappings
