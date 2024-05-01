import os

from rdkit import Chem

from gufe import SmallMoleculeComponent

root_path = os.path.dirname(__file__)
benzenes_sdf_path = f"{root_path}/benzenes_RHFE.sdf"

def get_benzene_compounds()->list[SmallMoleculeComponent]:
    """
    get the benzene test dataset parsed as SmallMoleculeComponents.

    Returns
    -------
    list[SmallMoleculeComponent]
        a list of the benzene compounds.

    """
    return [SmallMoleculeComponent.from_rdkit(rdm) for rdm in Chem.SDMolSupplier(benzenes_sdf_path, removeHs=False)]
