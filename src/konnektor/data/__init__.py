# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import os

from rdkit import Chem

from gufe import SmallMoleculeComponent

root_path = os.path.dirname(__file__)
benzenes_sdf_path = f"{root_path}/benzenes_RHFE.sdf"
hif2a_sdf_path = f"{root_path}/hif2a_ligands.sdf"
charged_ligands_path = f"{root_path}/charged_ligands.sdf"


def get_benzene_ligands() -> list[SmallMoleculeComponent]:
    """
    Get the benzene test dataset parsed as SmallMoleculeComponents.

    Returns
    -------
    list[SmallMoleculeComponent]
        a list of the benzene compounds.
    """
    return [
        SmallMoleculeComponent.from_rdkit(rdm)
        for rdm in Chem.SDMolSupplier(benzenes_sdf_path, removeHs=False)
    ]


def get_hif2a_ligands() -> list[SmallMoleculeComponent]:
    """
    Get the benzene test dataset parsed as SmallMoleculeComponents.

    Returns
    -------
    list[SmallMoleculeComponent]
        a list of the benzene compounds.
    """
    return [
        SmallMoleculeComponent.from_rdkit(rdm)
        for rdm in Chem.SDMolSupplier(hif2a_sdf_path, removeHs=False)
    ]


def get_charged_ligands() -> list[SmallMoleculeComponent]:
    """
    Get the benzene test dataset parsed as SmallMoleculeComponents.

    Returns
    -------
    list[SmallMoleculeComponent]
        a list of the benzene compounds.
    """
    return [
        SmallMoleculeComponent.from_rdkit(rdm)
        for rdm in Chem.SDMolSupplier(charged_ligands_path, removeHs=False)
    ]
