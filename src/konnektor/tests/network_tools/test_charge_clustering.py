import gufe
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from konnektor.network_tools.clustering import charge_clustering


@pytest.fixture
def formal_charge_mols():
    """Four simple molecules with +1, 0, 0, -1 formal charge."""

    def smiles_to_mol(smi):
        m = Chem.MolFromSmiles(smi)
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m)

        return m

    smiles = [
        "[NH3+]CCC(=O)O",
        "[NH2]CCC(=O)O",
        "[NH3+]CCC(=O)[O-]",
        "[NH2]CCC(=O)[O-]",
    ]

    mols = [smiles_to_mol(s) for s in smiles]

    return [gufe.SmallMoleculeComponent(m) for m in mols]


def test_charge_clustering(formal_charge_mols):
    clusterer = charge_clustering.ChargeClusterer()
    ret = clusterer.cluster_compounds(formal_charge_mols)

    assert isinstance(ret, dict)
    assert 0 in ret
    assert 1 in ret
    assert -1 in ret
    assert len(ret) == 3

    assert set(ret[1]) == {formal_charge_mols[0]}
    assert set(ret[0]) == {formal_charge_mols[1], formal_charge_mols[2]}
    assert set(ret[-1]) == {formal_charge_mols[-1]}
