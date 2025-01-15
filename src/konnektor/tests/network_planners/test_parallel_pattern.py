import pytest
import itertools
from gufe import ComponentMapping
from ...utils.toy_data import build_random_dataset
from ...network_planners.generators._parallel_mapping_pattern import (
    thread_mapping,
    _parallel_map_scoring,
)


def test_thread_mapping():
    components, mapper, scorer = build_random_dataset(n_compounds=20)
    component_pairs = list(itertools.combinations(components, 2))
    args = (1, component_pairs, [mapper], scorer)
    mappings = thread_mapping(args)

    assert len(mappings) == len(component_pairs)
    assert all(isinstance(m, ComponentMapping) for m in mappings)


@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
def test_parallel_map_scoring(with_progress, n_process):
    components, mapper, scorer = build_random_dataset(n_compounds=20)
    component_pairs = list(itertools.combinations(components, 2))

    mappings = _parallel_map_scoring(
        possible_edges=component_pairs,
        scorer=scorer,
        mappers=[mapper],
        n_processes=n_process,
        show_progress=with_progress,
    )

    assert len(mappings) == len(component_pairs)
    assert all(isinstance(m, ComponentMapping) for m in mappings)
