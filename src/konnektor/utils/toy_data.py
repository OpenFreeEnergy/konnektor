# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor


import numpy as np
from gufe import AtomMapper, AtomMapping, LigandAtomMapping, LigandNetwork, SmallMoleculeComponent
from rdkit import Chem


class emptyMapper(AtomMapper):
    def __init__(self):
        """
        Build a Mapper that only has use for dummy mappings. Generates empty mappings.
        """
        pass

    def suggest_mappings(self, molA, molB) -> AtomMapping:
        yield LigandAtomMapping(molA, molB, {})

    @classmethod
    def _defaults(cls):
        return super()._defaults()

    @classmethod
    def _from_dict(cls, d):
        s = cls()
        [setattr(s, k, v) for k, v in d.items()]
        return s

    def _to_dict(self):
        return vars(self)


def smiles_length_scorer(mapping: AtomMapping) -> float:
    int_A = len(mapping.componentA.smiles)
    int_B = len(mapping.componentB.smiles)

    return 1 / (int_A + int_B)


def build_random_dataset(n_compounds: int = 20, rand_seed: int = 13):
    """
    This function builds a random dataset of n_compounds artificial molecules.
    Additionally the generic scorer and mapper matching the compounds is returned.

    Parameters
    ----------
    n_compounds: int
        number of artificial molecules to build
    rand_seed: int
        random number seed.

    Returns
    -------
    (Iterable[SmallMoleculeComponent], AtomMapper, AtomMappingScorer)
        compounds, mapper, scorer
    """
    mapper = emptyMapper()
    scorer = smiles_length_scorer

    # generate random Molecules
    np.random.seed(rand_seed)
    smiles = [
        "".join(np.random.choice(["C", "O", "N", "S"], replace=True, size=(i % 10) + 1))
        for i in range(1, int(n_compounds) + 1)
    ]
    mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in smiles]
    [Chem.rdDistGeom.EmbedMolecule(m) for m in mols]

    compounds = [SmallMoleculeComponent(name=str(i), rdkit=m) for i, m in enumerate(mols)]

    return compounds, mapper, scorer


def build_random_mst_network(
    n_compounds=30, rand_seed=42, uni_score: bool = False
) -> LigandNetwork:
    """
    This function returns a randomized toy mst graph.

    Parameters
    ----------
    n_compounds:int
        number of artificial compounds
    rand_seed: int
        random seed number
    uni_score:
    whether to use always the score of 1 for an atom mapping

    Returns
    -------
    LigandNetwork
        the toy mst network
    """
    compounds, empty_mapper, scorer = build_random_dataset(
        n_compounds=n_compounds, rand_seed=rand_seed
    )

    if uni_score:
        scorer = lambda compound: 1

    from konnektor.network_planners import MinimalSpanningTreeNetworkGenerator

    planner = MinimalSpanningTreeNetworkGenerator(mappers=empty_mapper, scorer=scorer)

    ligand_network = planner(compounds)
    return ligand_network


def build_n_random_mst_network(
    n_compounds=30,
    rand_seed=42,
    sub_networks: int = 2,
    overlap: int = 1,
    uni_score: bool = False,
) -> tuple[LigandNetwork, LigandNetwork]:
    """
    This function returns a randomized toy mst graph.

    Parameters
    ----------
    n_compounds:int
        number of artificial compounds
    rand_seed: int
        random seed number
    uni_score:
    whether to use always the score of 1 for an atom mapping

    Returns
    -------
    LigandNetwork
        the toy mst network
    """
    compounds, empty_mapper, scorer = build_random_dataset(
        n_compounds=n_compounds, rand_seed=rand_seed
    )

    if uni_score:
        scorer = lambda compound: 1

    from konnektor.network_planners import MinimalSpanningTreeNetworkGenerator

    planner = MinimalSpanningTreeNetworkGenerator(mappers=empty_mapper, scorer=scorer)

    networks = []
    step = n_compounds // sub_networks
    for i in range(sub_networks):
        if i == sub_networks - 1:
            sub_components = compounds[i * step - overlap :]
        elif i > 0:
            sub_components = compounds[i * step - overlap : (i + 1) * step]
        else:
            sub_components = compounds[i * step : (i + 1) * step]
        networks.append(planner(sub_components))

    return networks


def build_random_fully_connected_network(
    n_compounds=30, rand_seed=42, uni_score: bool = False
) -> LigandNetwork:
    """
    This function returns a randomized toy fully connected graph.

    Parameters
    ----------
    n_compounds:int
        number of artificial compounds
    rand_seed: int
        random seed number
    uni_score:
        whether to use always the score of 1 for an atom mapping

    Returns
    -------
    LigandNetwork
        the toy fully connected network
    """
    compounds, empty_mapper, scorer = build_random_dataset(
        n_compounds=n_compounds, rand_seed=rand_seed
    )

    if uni_score:
        scorer = lambda compound: 1

    from konnektor.network_planners import MaximalNetworkGenerator

    planner = MaximalNetworkGenerator(mappers=empty_mapper, scorer=scorer)

    ligand_network = planner(compounds)
    return ligand_network
