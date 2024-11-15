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


class DummyMapper(AtomMapper):
    def __init__(self, id:int):
        """
        Build a generic Mapper, that only has use for dummy mappings.
        Generates mappings that increment by 1 for testing order-dependent behavior.

        id : use to identify an instance of this class in testing
        """
        self.id = id

    def suggest_mappings(self, molA, molB) -> AtomMapping:
        for i in range(3):
            yield LigandAtomMapping(molA, molB, {int(self.id):i})

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

def test_generate_maximal_network_missing_scorer(toluene_vs_others):
    """If no scorer is provided, the first mapping of the last mapper should be used."""

    toluene, others = toluene_vs_others
    components = others+[toluene]

    id_a = 0
    id_b = 1
    planner = MaximalNetworkGenerator(
        mappers= [DummyMapper(id=id_a), DummyMapper(id=id_b)],
        scorer=None,
        progress=False,
        n_processes=1,
    )

    network = planner.generate_ligand_network(components)

    assert [e.componentA_to_componentB for e in network.edges] == len(network.edges)*[{id_b:0}]
