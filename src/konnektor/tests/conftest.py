# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe
from collections.abc import Iterable
from typing import NamedTuple

import pytest
from gufe import LigandAtomMapping, LigandNetwork, SmallMoleculeComponent

from .network_planners.conf import mol_from_smiles


@pytest.fixture
def mols():
    mol1 = SmallMoleculeComponent(mol_from_smiles("CCO"))
    mol2 = SmallMoleculeComponent(mol_from_smiles("CC"))
    mol3 = SmallMoleculeComponent(mol_from_smiles("CO"))
    return mol1, mol2, mol3


@pytest.fixture
def std_edges(mols):
    mol1, mol2, mol3 = mols
    edge12 = LigandAtomMapping(mol1, mol2, {0: 0, 1: 1})
    edge23 = LigandAtomMapping(mol2, mol3, {0: 0})
    edge13 = LigandAtomMapping(mol1, mol3, {0: 0, 2: 1})
    return edge12, edge23, edge13


class _NetworkTestContainer(NamedTuple):
    """Container to facilitate network testing"""

    network: LigandNetwork
    nodes: Iterable[SmallMoleculeComponent]
    edges: Iterable[LigandAtomMapping]
    n_nodes: int
    n_edges: int


@pytest.fixture
def simple_network(mols, std_edges):
    """Network with no edges duplicated and all nodes in edges"""
    network = LigandNetwork(std_edges)
    return _NetworkTestContainer(
        network=network,
        nodes=mols,
        edges=std_edges,
        n_nodes=3,
        n_edges=3,
    )
