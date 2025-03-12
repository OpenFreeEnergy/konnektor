import importlib

import pytest
from gufe import LigandAtomMapping, LigandNetwork, SmallMoleculeComponent
from rdkit import Chem


@pytest.fixture(scope="session")
def toluene_vs_others(atom_mapping_basic_test_files):
    central_ligand_name = "toluene"
    others = [v for (k, v) in atom_mapping_basic_test_files.items() if k != central_ligand_name]
    toluene = atom_mapping_basic_test_files[central_ligand_name]
    return toluene, others


@pytest.fixture(scope="session")
def atom_mapping_basic_test_files():
    # a dict of {filenames.strip(mol2): SmallMoleculeComponent} for a simple
    # set of ligands
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
def ligand_network_ab(atom_mapping_basic_test_files):
    mappingsA = []
    mappingsB = []

    for _, mA in atom_mapping_basic_test_files.items():
        for _, mB in atom_mapping_basic_test_files.items():
            m = LigandAtomMapping(mA, mB, {}).with_annotations({"score": 1})
            if mA == mB:
                continue
            elif mA.name.startswith("2") and mB.name.startswith("2"):
                mappingsA.append(m)
            elif not (mA.name.startswith("2") or mB.name.startswith("2")):
                mappingsB.append(m)

    lna = LigandNetwork(edges=mappingsA)
    lnb = LigandNetwork(edges=mappingsB)
    return (lna, lnb)
