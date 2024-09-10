# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from rdkit import Chem

try:
    # R-group enumeration intermediate generator
    from rgroupinterm import rgroupenumeration
    from rgroupinterm import pruners
except:
    pass

from gufe import SmallMoleculeComponent

from ._abstract_intermediator import Intermediator
from ...utils.optional_import import requires_package


@requires_package("rgroupinterm")
class ImergeIntermediator(Intermediator):
    def __init__(
        self,
        enumerate_kekule: bool = False,
        permutate: bool = False,
        insert_small: bool = False,
    ):
        """

        Parameters
        ----------
        enumerate_kekule: bool, optional
            (default: False)
        permutate: bool, optional
            (default: False)
        insert_small: bool, optional
            (default: False)
        """
        self.enumerate_kekule = enumerate_kekule
        self.permutate = permutate
        self.insert_small = insert_small

    def generate_intermediate(
        self, molA: SmallMoleculeComponent, molB: SmallMoleculeComponent
    ) -> SmallMoleculeComponent:
        """
        Uses the rgroupinterm package to generate a intermediate between molA and molB.

        Parameters
        ----------
        molA: SmallMoleculeComponent
        molB: SmallMoleculeComponent

        Returns
        -------
        Iterator[SmallMoleculeComponent]
            returns the small molecule intermediate between molA and molB.

        """

        rdmolA = Chem.MolFromSmiles(Chem.MolToSmiles(molA.to_rdkit()))
        rdmolB = Chem.MolFromSmiles(Chem.MolToSmiles(molB.to_rdkit()))

        #  get rgroups
        generator = rgroupenumeration.EnumRGroups(
            self.enumerate_kekule, self.permutate, self.insert_small
        )
        df_interm, _ = generator.generate_intermediates([rdmolA, rdmolB])

        # prune groups to get intermediates
        pruner = pruners.BasePruner(
            [
                pruners.TanimotoScorer(
                    transformer=pruners.HarmonicMeanTransformer(exponent=4)
                )
            ],
            topn=1,
        )

        df_interm["Parent_1"] = rdmolA
        df_interm["Parent_2"] = rdmolB
        df_interm["Pair"] = 0
        pruned_df = pruner(df_interm)

        # intermediate generation
        rd_mol_intermediate = pruned_df["Intermediate"].values[0]
        rd_mol_intermediate = Chem.AddHs(rd_mol_intermediate)
        Chem.AllChem.EmbedMolecule(rd_mol_intermediate)

        mol_intermediate = SmallMoleculeComponent.from_rdkit(
            rdkit=rd_mol_intermediate,
            name=molA.name + "_" + molB.name + "_intermediate",
        )

        yield mol_intermediate
