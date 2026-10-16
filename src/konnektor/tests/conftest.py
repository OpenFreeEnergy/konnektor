# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe
import importlib
from collections.abc import Iterable
from typing import NamedTuple

import pytest
from gufe import LigandAtomMapping, LigandNetwork, SmallMoleculeComponent
from rdkit import Chem

from .network_planners.conf import mol_from_smiles


@pytest.fixture(scope="session")
def atom_mapping_basic_test_files():
    """a dict of {filenames.strip(mol2): SmallMoleculeComponent} for a simple set of ligands"""

    files = {}
    for f in [
        "1,3,7-trimethylnaphthalene",
        "1-butyl-4-methylbenzene",
        "2,6-dimethylnaphthalene",
        "2-methyl-6-propylnaphthalene",
        "2-methylnaphthalene",
        "2-naftanol",
        "methylcyclohexane",
        "toluene",
    ]:
        with importlib.resources.path("konnektor.tests.data", f + ".mol2") as fn:
            mol = Chem.MolFromMol2File(str(fn), removeHs=False)
            files[f] = SmallMoleculeComponent(mol, name=f)

    return files


@pytest.fixture(scope="session")
def toluene_vs_others(atom_mapping_basic_test_files):
    central_ligand_name = "toluene"
    others = [v for (k, v) in atom_mapping_basic_test_files.items() if k != central_ligand_name]
    toluene = atom_mapping_basic_test_files[central_ligand_name]
    return toluene, others


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
