# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

from konnektor.network_planners import MaximalNetworkGenerator
from konnektor.tests.network_planners.conf import (
    atom_mapping_basic_test_files,
    toluene_vs_others,
    genScorer,
    GenAtomMapper,
)
from gufe import LigandAtomMapping, AtomMapper, AtomMapping
from konnektor.utils.toy_data import build_random_dataset


class CounterTestMapper(AtomMapper):
    def __init__(self):
        """
        Build a generic Mapper, that only has use for dummy mappings.
        Generates mappings that increment by 1 for testing order-dependent behavior.
        """
        pass

    def suggest_mappings(self, molA, molB) -> AtomMapping:
        yield LigandAtomMapping(molA, molB, {0:1})

    @classmethod
    def _defaults(cls):
        return super()._defaults()

    @classmethod
    def _from_dict(cls, d):
        s = cls()
        [setattr(s, k, v) for k, v in d.items()]
        return s

    def _to_dict(self):
        return vars(self)


@pytest.mark.parametrize("n_process", [1, 2])
@pytest.mark.parametrize("with_progress", [True, False])
@pytest.mark.parametrize("with_scorer", [True, False])
def test_generate_maximal_network(
    toluene_vs_others, with_progress, with_scorer, n_process
):
    toluene, others = toluene_vs_others

    mapper = GenAtomMapper()

    scorer = genScorer if with_scorer else None

    planner = MaximalNetworkGenerator(
        mappers=mapper, scorer=scorer, progress=with_progress, n_processes=n_process
    )
    network = planner.generate_ligand_network(others + [toluene])

    assert len(network.nodes) == len(others) + 1

    edge_count = len(others) * (len(others) + 1) / 2
    assert len(network.edges) == edge_count

    if scorer:
        for edge in network.edges:
            score = edge.annotations["score"]
            assert score == 1.0 / len(edge.componentA_to_componentB)
    else:
        for edge in network.edges:
            assert "score" not in edge.annotations

def test_generate_maximal_network_missing_scorer():
    """If no scorer is provided, the last mapper tried should be used."""

    components, empty_mapper, _ = build_random_dataset(n_compounds=2)
    counter_mapper = CounterTestMapper()
    planner = MaximalNetworkGenerator(
        mappers= [empty_mapper, counter_mapper],
        scorer=None,
        progress=False,
        n_processes=1,
    )
    network = planner.generate_ligand_network(components)
    assert [e for e in network.edges][-1].componentA_to_componentB == {0:1}
    [(e.componentA.name, e.componentB.name) for e in network.edges]