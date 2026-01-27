# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import functools
import multiprocessing as mult
import warnings
from collections.abc import Callable

from gufe import AtomMapper, AtomMapping, SmallMoleculeComponent
from tqdm.auto import tqdm


def _determine_best_mapping(
    component_pair: tuple[SmallMoleculeComponent, SmallMoleculeComponent],
    mappers: list[AtomMapper],
    scorer: Callable | None,
) -> AtomMapping:
    """
    Helper function to generate all possible mappings for ``component_pair`` using the ``mappers``.
    If a ``scorer`` is provided, score all mappings and return the mapping with the best score.
    If no ``scorer`` is provided, return the first mapping from the first mapper.


    Parameters
    ----------
    component_pair : tuple[SmallMoleculeComponent, SmallMoleculeComponent]
        The two molecules for which the best mapping will be determined.
    mappers : AtomMapper | list[AtomMapper]
        The mapper(s) to use to generate possible mappings between the molecules in the ``component_pair``.
    scorer : Optional[Callable]
        The mapping scorer to use, in the form of a ``Callable`` that takes in an ``AtomMapping`` and returns a float in [0,1].

    Returns
    -------
    AtomMapping
        The best mapping found for the component pair.
    """
    best_score = 0.0
    best_mapping = None
    molA = component_pair[0]
    molB = component_pair[1]

    for mapper in mappers:
        try:
            mapping_generator = mapper.suggest_mappings(molA, molB)
        except:  # TODO: I don't like this bare except
            continue

        if scorer:
            tmp_mappings = [
                mapping.with_annotations({"score": scorer(mapping)})
                for mapping in mapping_generator
            ]

            if len(tmp_mappings) > 0:
                # TODO: where should we enforce that this score is in [0,1]?
                tmp_best_mapping = max(tmp_mappings, key=lambda m: m.annotations["score"])
                # TODO: we still need a more explicit tie-breaking scheme
                if tmp_best_mapping.annotations["score"] > best_score or best_mapping is None:
                    best_score = tmp_best_mapping.annotations["score"]
                    best_mapping = tmp_best_mapping
        else:
            try:
                best_mapping = next(mapping_generator)
                if len(mappers) > 1:
                    warnings.warn(
                        "Multiple mappers were provided, but no scorer. "
                        f"Only the first valid mapper will be used: {mapper}"
                    )
                break
            except:  # TODO: fix this bare except, or remove it (first mapper vs. first valid mapper)
                continue

    return best_mapping


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
        best_mapping = _determine_best_mapping(
            component_pair=component_pair, mappers=mappers, scorer=scorer
        )
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
    Parallelize mapping and scoring of a given list of molecule pairs.

    Parameters
    ----------
    possible_edges: tuple[SmallMoleculeComponent, SmallMoleculeComponent]
        two  molecules to be mapped.
    scorer: callable
        scoring the mappings
    mappers: AtomMapper
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


def _serial_map_scoring(
    possible_edges: list[tuple[SmallMoleculeComponent, SmallMoleculeComponent]],
    scorer: Callable[[AtomMapping], float],
    mappers: list[AtomMapper],
    show_progress: bool = True,
) -> list[AtomMapping]:
    """_summary_

    Parameters
    ----------
    possible_edges: tuple[SmallMoleculeComponent, SmallMoleculeComponent]
        two  molecules to be mapped.
    scorer: callable
        scoring the mappings
    mappers: AtomMapper
        atom mapper for the mappings
    show_progress: bool
        show a tqdm progressbar.

    Returns
    -------
    list[AtomMapping]:
        return a list of scored atom mappings
    """
    if show_progress is True:
        progress = functools.partial(tqdm, total=len(possible_edges), delay=1.5, desc="Mapping")
    else:
        progress = lambda x: x

    mappings = []
    for component_pair in progress(possible_edges):
        best_mapping = _determine_best_mapping(
            component_pair=component_pair, mappers=mappers, scorer=scorer
        )

        if best_mapping is not None:
            mappings.append(best_mapping)

    return mappings


def _score_mappings(n_processes: int, **kwargs) -> list[AtomMapping]:
    if n_processes > 1:
        scored_mappings = _parallel_map_scoring(
            n_processes=n_processes,
            **kwargs,
        )
    else:  # serial variant
        scored_mappings = _serial_map_scoring(**kwargs)

    return scored_mappings
