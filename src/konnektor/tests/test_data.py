# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from gufe import SmallMoleculeComponent
from konnektor.data import get_benzene_ligands, get_hif2a_ligands, get_charged_ligands


def test_get_benzenes():
    compounds = get_benzene_ligands()

    assert len(compounds) == 16
    assert all(isinstance(c, SmallMoleculeComponent) for c in compounds)


def test_get_hif2a():
    compounds = get_hif2a_ligands()

    assert len(compounds) == 37
    assert all(isinstance(c, SmallMoleculeComponent) for c in compounds)


def test_get_hif2a():
    compounds = get_charged_ligands()

    assert len(compounds) == 44
    assert all(isinstance(c, SmallMoleculeComponent) for c in compounds)
