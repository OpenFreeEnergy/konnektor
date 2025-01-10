# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

from rdkit import Chem

from gufe import SmallMoleculeComponent
from ...network_tools import ImergeIntermediator


@pytest.fixture
def example_mols():
    rdmolA = Chem.AddHs(
        Chem.MolFromSmiles("CC(C)(C)C1=CC=C2N[C@@H](C3=CC=C(O)C=C3)[C@@H]3CCCO[C@@H]3C2=C1")
    )
    Chem.rdDistGeom.EmbedMolecule(rdmolA)
    molA = SmallMoleculeComponent.from_rdkit(rdmolA, name="molA")
    rdmolB = Chem.AddHs(
        Chem.MolFromSmiles("[NH3+]C[C@H]1CC[C@@H]2[C@H](O1)C1=CC(C(F)(F)F)=CC=C1N[C@H]2C1=CC=CC=C1")
    )
    Chem.rdDistGeom.EmbedMolecule(rdmolB)
    molB = SmallMoleculeComponent.from_rdkit(rdmolB, name="molB")
    return molA, molB


@pytest.importorskip("ImergeIntermediator")
def test_imerge_imediator_gen_intermediate(example_mols):
    molA, molB = example_mols

    expected_rdmol = Chem.AddHs(
        Chem.MolFromSmiles("CC(C)(C)c1ccc2c(c1)[C@H]1O[C@@H](C[NH3+])CC[C@H]1[C@H](c1ccc(O)cc1)N2")
    )
    Chem.rdDistGeom.EmbedMolecule(expected_rdmol)
    expected_mol = SmallMoleculeComponent.from_rdkit(expected_rdmol, name="intermediate")

    intermediator = ImergeIntermediator()
    intermediate_mol = next(intermediator(molA, molB))

    assert expected_mol == intermediate_mol
