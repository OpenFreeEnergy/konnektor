import importlib

import pytest
from gufe import SmallMoleculeComponent
from rdkit import Chem


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
