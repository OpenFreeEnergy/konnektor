# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import inspect

from gufe import AtomMapper, LigandAtomMapping, SmallMoleculeComponent
from rdkit import Chem
from rdkit.Chem import AllChem


def mol_from_smiles(smiles: str) -> Chem.Mol:
    m = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(m)

    return m


class DummyAtomMapper(AtomMapper):
    @classmethod
    def _from_dict(cls, d: dict):
        """Deserialize from dict representation"""
        if any(k not in cls._defaults() for k in d):
            keys = list(filter(lambda k: k in cls._defaults(), d.keys()))
            raise ValueError(f"I don't know about all the keys here: {keys}")
        return cls(**d)

    def _to_dict(self) -> dict:
        d = {}
        for key in self._defaults():
            if hasattr(self, key):
                d[key] = getattr(self, key)
        return d

    @classmethod
    def _defaults(cls):
        """This method should be overridden to provide the dict of defaults
        appropriate for the `GufeTokenizable` subclass.
        """
        sig = inspect.signature(cls.__init__)

        defaults = {
            param.name: param.default
            for param in sig.parameters.values()
            if param.default is not inspect.Parameter.empty
        }

        return defaults


class GenAtomMapper(DummyAtomMapper):
    """Generic mapper that makes a single valid mapping."""
    def suggest_mappings(
        self, componentA: SmallMoleculeComponent, componentB: SmallMoleculeComponent
    ):
        atomsA = range(componentA.to_rdkit().GetNumAtoms())
        atomsB = range(componentB.to_rdkit().GetNumAtoms())

        yield LigandAtomMapping(
            componentA,
            componentB,
            componentA_to_componentB={k: v for k, v in zip(atomsA, atomsB)},
        )


class MultiAtomMapper(DummyAtomMapper):
    """Returns two mapping suggestions, the second mapping is just reversed"""
    def suggest_mappings(
        self, componentA: SmallMoleculeComponent, componentB: SmallMoleculeComponent
    ):
        atomsA = range(componentA.to_rdkit().GetNumAtoms())
        atomsB = list(range(componentB.to_rdkit().GetNumAtoms()))

        # first map to atomsB forwards, then map in the reverse order
        for dir in [1, -1]:
            yield LigandAtomMapping(
                componentA,
                componentB,
                componentA_to_componentB={k: v for k, v in zip(atomsA, atomsB[::dir])},
            )


class BadMapper(DummyAtomMapper):
    def suggest_mappings(
        self, componentA: SmallMoleculeComponent, componentB: SmallMoleculeComponent
    ):
        yield LigandAtomMapping(componentA, componentB, componentA_to_componentB={0: 0})


class SuperBadMapper(DummyAtomMapper):
    def suggest_mappings(
        self, componentA: SmallMoleculeComponent, componentB: SmallMoleculeComponent
    ):
        yield LigandAtomMapping(componentA, componentB, componentA_to_componentB={})


class ErrorMapper(DummyAtomMapper):
    def suggest_mappings(
        self, componentA: SmallMoleculeComponent, componentB: SmallMoleculeComponent
    ):
        raise ValueError("No mapping found for")


def genScorer(mapping):
    return 1.0 / len(mapping.componentA_to_componentB)
