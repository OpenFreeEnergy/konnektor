# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
from gufe import LigandNetwork

from konnektor.network_planners.concatenators._abstract_network_concatenator import (
    NetworkConcatenator,
)
from konnektor.utils.toy_data import (
    EmptyMapper,
    RandomScorer,
    build_random_dataset,
)


class DummyConcatenator(NetworkConcatenator):
    def __init__(self, result=None):
        super().__init__(
            mappers=EmptyMapper(),
            scorer=RandomScorer(n=2),
            network_generator=None,
        )
        self.result = result
        self.received_exclude = None

    def _concatenate_networks(
        self,
        ligand_networks,
        exclude,
    ):
        self.received_exclude = exclude

        if self.result is not None:
            return self.result

        return ligand_networks[0]


@pytest.fixture
def two_ligands():
    components, _, _ = build_random_dataset(
        n_compounds=2,
        rand_seed=1,
    )
    return components


def test_concatenate_networks_requires_input():
    concatenator = DummyConcatenator()

    with pytest.raises(
        ValueError,
        match="At least one LigandNetwork is required",
    ):
        concatenator.concatenate_networks([])


def test_concatenate_networks_rejects_disconnected_input(two_ligands):
    ligand_a, ligand_b = two_ligands
    network = LigandNetwork(nodes=[ligand_a, ligand_b], edges=[])
    concatenator = DummyConcatenator()

    with pytest.raises(
        RuntimeError,
        match="input networks are disconnected",
    ):
        concatenator.concatenate_networks([network])


def test_concatenate_networks_returns_single_network_unchanged(two_ligands):
    ligand_a, _ = two_ligands
    network = LigandNetwork(nodes=[ligand_a], edges=[])
    concatenator = DummyConcatenator()

    result = concatenator.concatenate_networks([network])

    assert result is network


def test_concatenate_networks_normalizes_excluded_edges(two_ligands):
    ligand_a, ligand_b = two_ligands

    network_a = LigandNetwork(nodes=[ligand_a], edges=[])
    network_b = LigandNetwork(nodes=[ligand_b], edges=[])
    concatenator = DummyConcatenator()
    mapping = next(EmptyMapper().suggest_mappings(ligand_a, ligand_b))

    concatenator.concatenate_networks(
        [network_a, network_b],
        exclude_edges=[mapping],
    )

    assert concatenator.received_exclude == {frozenset((mapping.componentA, mapping.componentB))}


def test_concatenate_networks_rejects_disconnected_result(two_ligands):
    ligand_a, ligand_b = two_ligands

    network_a = LigandNetwork(nodes=[ligand_a], edges=[])
    network_b = LigandNetwork(nodes=[ligand_b], edges=[])
    disconnected_result = LigandNetwork(nodes=[ligand_a, ligand_b], edges=[])
    concatenator = DummyConcatenator(result=disconnected_result)

    with pytest.raises(RuntimeError, match="Could not build a connected network"):
        concatenator.concatenate_networks([network_a, network_b])
