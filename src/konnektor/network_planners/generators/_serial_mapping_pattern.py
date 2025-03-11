import functools
from collections.abc import Callable

from gufe import AtomMapper, AtomMapping, SmallMoleculeComponent
from tqdm.auto import tqdm


def _serial_map_scoring(
    possible_edges: list[tuple[SmallMoleculeComponent, SmallMoleculeComponent]],
    scorer: Callable[[AtomMapping], float],
    mappers: list[AtomMapper],
    edges_to_score: int,
    show_progress: bool = True,
):
    if show_progress is True:
        progress = functools.partial(tqdm, total=edges_to_score, delay=1.5, desc="Mapping")
    else:
        progress = lambda x: x

    mappings = []
    for component_pair in progress(possible_edges):
        best_score = 0.0
        best_mapping = None
        molA = component_pair[0]
        molB = component_pair[1]

        for mapper in mappers:
            try:
                mapping_generator = mapper.suggest_mappings(molA, molB)
            except:
                continue

            if scorer:
                tmp_mappings = [
                    mapping.with_annotations({"score": scorer(mapping)})
                    for mapping in mapping_generator
                ]

                if len(tmp_mappings) > 0:
                    tmp_best_mapping = min(tmp_mappings, key=lambda m: m.annotations["score"])

                    if tmp_best_mapping.annotations["score"] < best_score or best_mapping is None:
                        best_score = tmp_best_mapping.annotations["score"]
                        best_mapping = tmp_best_mapping
            else:
                try:
                    best_mapping = next(mapping_generator)
                except:
                    print("warning")
                    continue

        if best_mapping is not None:
            mappings.append(best_mapping)

    return mappings
