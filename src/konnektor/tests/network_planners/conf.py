# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import importlib
import pytest
from importlib import resources
from rdkit import Chem
from rdkit.Chem import AllChem

import konnektor
from konnektor import network_planners

from gufe import AtomMapper
from gufe import SmallMoleculeComponent, LigandAtomMapping
from gufe import LigandNetwork, LigandAtomMapping



def mol_from_smiles(smiles: str) -> Chem.Mol:
    m = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(m)

    return m

@pytest.fixture(scope='session')
def atom_mapping_basic_test_files():
    # a dict of {filenames.strip(mol2): SmallMoleculeComponent} for a simple
    # set of ligands
    files = {}
    for f in [
        '1,3,7-trimethylnaphthalene',
        '1-butyl-4-methylbenzene',
        '2,6-dimethylnaphthalene',
        '2-methyl-6-propylnaphthalene',
        '2-methylnaphthalene',
        '2-naftanol',
        'methylcyclohexane',
        'toluene']:
        with importlib.resources.path('konnektor.tests.data',
                                      f + '.mol2') as fn:
            mol = Chem.MolFromMol2File(str(fn), removeHs=False)
            files[f] = SmallMoleculeComponent(mol, name=f)

    return files

class GenAtomMapper(AtomMapper):
    def suggest_mappings(self, componentA:SmallMoleculeComponent,
                         componentB:SmallMoleculeComponent):
        atomsA = range(componentA.to_rdkit().GetNumAtoms())
        atomsB = range(componentB.to_rdkit().GetNumAtoms())

        yield LigandAtomMapping(componentA, componentB,
                                componentA_to_componentB={k:v for k,v in zip(
                                    atomsA, atomsB)})

class BadMapper(AtomMapper):
    def suggest_mappings(self, componentA:SmallMoleculeComponent,
                         componentB:SmallMoleculeComponent):
        yield LigandAtomMapping(componentA, componentB,
                                componentA_to_componentB={0:0})

class SuperBadMapper(AtomMapper):
    def suggest_mappings(self, componentA:SmallMoleculeComponent,
                         componentB:SmallMoleculeComponent):
        yield LigandAtomMapping(componentA, componentB,
                                componentA_to_componentB={})


class ErrorMapper(AtomMapper):
    def suggest_mappings(self, componentA:SmallMoleculeComponent,
                         componentB:SmallMoleculeComponent):
        #raise StopIteration('No mapping found for')# Check for good solution
        # here
        raise ValueError('No mapping found for')

def genScorer(mapping):
    return 1.0 / len(mapping.componentA_to_componentB)


@pytest.fixture(scope='session')
def toluene_vs_others(atom_mapping_basic_test_files):
    central_ligand_name = 'toluene'
    others = [v for (k, v) in atom_mapping_basic_test_files.items()
              if k != central_ligand_name]
    toluene = atom_mapping_basic_test_files[central_ligand_name]
    return toluene, others



@pytest.fixture(scope='session')
def ligand_network_ab(atom_mapping_basic_test_files):
        mappingsA = []
        mappingsB = []

        for _,mA in atom_mapping_basic_test_files.items():
            for _,mB in atom_mapping_basic_test_files.items():
                m= LigandAtomMapping(mA, mB, {}
                                     ).with_annotations({'score': 1})
                if mA==mB:
                    continue
                elif mA.name.startswith("2") and mB.name.startswith("2"):
                    mappingsA.append(m)
                elif not (mA.name.startswith("2") or mB.name.startswith("2")):
                    mappingsB.append(m)

        lna = LigandNetwork(edges=mappingsA)
        lnb = LigandNetwork(edges=mappingsB)
        return (lna, lnb)