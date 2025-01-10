# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from typing import Tuple

import numpy as np
from gufe import LigandAtomMapping, AtomMapper, AtomMapping
from gufe import SmallMoleculeComponent, LigandNetwork
from rdkit import Chem


class genMapper(AtomMapper):
    def __init__(self):
        """
        Build a generic Mapper, that only has use for dummy mappings. Generates empty mappings
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


class genScorer:  # (AtomMappingScorer):
    def __init__(self, n_scores: int, rand_seed: int = None):
        """
        Builds a scorer that contains a predefined sequence of scores, n_scores long and each score is initially randomly uniformly picked between 1 and 0.
        The scorer repeats the score sequence after n_scoresth time, calling the scorer obj.
        The use of this class is currently envisioned for toydata and testing.

        Parameters
        ----------
        n_scores: int
            number of scores to build
        rand_seed: int
            random number seed for the random scores.
        """
        np.random.seed(rand_seed)

        self.vals = np.random.uniform(size=n_scores)

        self.vals = self.vals.astype(float)
        self.n_scores = n_scores
        self.i = 0

    def __call__(self, mapping):
        # todo: remove once subclassed from gufe
        return self.get_score(mapping)

    def get_score(self, mapping: AtomMapping) -> float:
        """
        return the score, at position self.i

        Parameters
        ----------
        mapping: AtomMapping
            the score will not be depending on the mapping! this mimicks only classical scorer use.

        Returns
        -------
        float
            score to be returned.

        """
        v = self.vals[self.i]
        self.i = (self.i + 1) % self.n_scores
        return v


def build_random_dataset(n_compounds: int = 20, rand_seed: int = None):
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
    gen_mapper = genMapper()
    gen_scorer = genScorer(n_scores=n_compounds, rand_seed=rand_seed)

    # generate random Molecules
    np.random.seed(rand_seed)
    smiles = [
        "".join(np.random.choice(["C", "O", "N", "S"], replace=True, size=(i % 30) + 1))
        for i in range(1, int(n_compounds) + 1)
    ]
    mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in smiles]
    [Chem.rdDistGeom.EmbedMolecule(m) for m in mols]
    compounds = [SmallMoleculeComponent(name=str(i), rdkit=m) for i, m in enumerate(mols)]

    return compounds, gen_mapper, gen_scorer


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
    compounds, genMapper, genScorer = build_random_dataset(
        n_compounds=n_compounds, rand_seed=rand_seed
    )

    if uni_score:
        genScorer.get_score = lambda compound: 1

    from konnektor.network_planners import MinimalSpanningTreeNetworkGenerator

    planner = MinimalSpanningTreeNetworkGenerator(mappers=genMapper, scorer=genScorer)

    ligand_network = planner(compounds)
    return ligand_network


def build_n_random_mst_network(
    n_compounds=30,
    rand_seed=42,
    sub_networks: int = 2,
    overlap: int = 1,
    uni_score: bool = False,
) -> Tuple[LigandNetwork, LigandNetwork]:
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
    compounds, genMapper, genScorer = build_random_dataset(
        n_compounds=n_compounds, rand_seed=rand_seed
    )

    if uni_score:
        genScorer.get_score = lambda compound: 1

    from konnektor.network_planners import MinimalSpanningTreeNetworkGenerator

    planner = MinimalSpanningTreeNetworkGenerator(mappers=genMapper, scorer=genScorer)

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
    compounds, genMapper, genScorer = build_random_dataset(
        n_compounds=n_compounds, rand_seed=rand_seed
    )

    if uni_score:
        genScorer.get_score = lambda compound: 1

    from konnektor.network_planners import MaximalNetworkGenerator

    planner = MaximalNetworkGenerator(mappers=genMapper, scorer=genScorer)

    ligand_network = planner(compounds)
    return ligand_network
