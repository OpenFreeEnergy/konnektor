import pytest
from gufe import LigandAtomMapping, LigandNetwork


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
