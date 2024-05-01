
import pytest
from gufe import SmallMoleculeComponent
from konnektor.data import get_benzene_compounds

def test_get_benzenes():
    compounds = get_benzene_compounds()

    assert len(compounds) == 16
    assert all(isinstance(c, SmallMoleculeComponent) for c in compounds)
