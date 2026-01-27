import itertools

import pytest
from gufe import ComponentMapping

from konnektor.network_planners._map_scoring import (
    _parallel_map_scoring,
    _serial_map_scoring,
    thread_mapping,
)

from ...utils.toy_data import build_random_dataset


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


def test_parallel_serial_equality():
    components, mapper, scorer = build_random_dataset(n_compounds=20)
    component_pairs = list(itertools.combinations(components, 2))

    mappings_parallel = _parallel_map_scoring(
        possible_edges=component_pairs,
        scorer=scorer,
        mappers=[mapper],
        n_processes=1,
        show_progress=False,
    )
    mappings_serial = _serial_map_scoring(
        possible_edges=component_pairs,
        scorer=scorer,
        mappers=[mapper],
        n_edges_to_score=len(component_pairs),
        show_progress=False,
    )
    assert len(mappings_parallel) == len(mappings_serial)
    comps_serial = {(m.componentA, m.componentB) for m in mappings_serial}
    comps_parallel = {(m.componentA, m.componentB) for m in mappings_parallel}
    assert comps_serial == comps_parallel
    assert set(mappings_parallel) == set(mappings_serial)
    assert len(mappings_parallel) == len(component_pairs)
    assert all(isinstance(m, ComponentMapping) for m in mappings_parallel)
