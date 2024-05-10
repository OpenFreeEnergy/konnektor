import numpy as np
import pytest
from gufe import AtomMapper, AtomMapping, SmallMoleculeComponent, LigandNetwork
from rdkit import Chem

from konnektor.utils.toy_data import genMapper, genScorer, build_random_dataset, \
    build_random_mst_network, \
    build_random_fully_connected_network, build_n_random_mst_network


def test_genMapper():
    rdmolA = Chem.MolFromSmiles('c1ccccc1')
    rdmolA = Chem.AddHs(rdmolA)
    Chem.rdDistGeom.EmbedMolecule(rdmolA)
    molA = SmallMoleculeComponent.from_rdkit(rdmolA)

    rdmolB = Chem.MolFromSmiles('Cc1ccccc1')
    rdmolB = Chem.AddHs(rdmolB)
    Chem.rdDistGeom.EmbedMolecule(rdmolB)
    molB = SmallMoleculeComponent.from_rdkit(rdmolB)

    mapper = genMapper()
    mapping = next(mapper.suggest_mappings(molA, molB))

    assert isinstance(mapper, AtomMapper)
    assert isinstance(mapping, AtomMapping)
    assert hasattr(mapping, "componentA_to_componentB")


def test_genScorer():
    n_scores = 10
    scorer = genScorer(n_scores=n_scores, rand_seed=42)

    # assert isinstance(scorer, AtomMappingScorer)
    assert len(scorer.vals) == n_scores
    assert scorer.i == 0

    # get some input data
    rdmolA = Chem.MolFromSmiles('c1ccccc1')
    rdmolA = Chem.AddHs(rdmolA)
    Chem.rdDistGeom.EmbedMolecule(rdmolA)
    molA = SmallMoleculeComponent.from_rdkit(rdmolA)

    rdmolB = Chem.MolFromSmiles('Cc1ccccc1')
    rdmolB = Chem.AddHs(rdmolB)
    Chem.rdDistGeom.EmbedMolecule(rdmolB)
    molB = SmallMoleculeComponent.from_rdkit(rdmolB)

    mapper = genMapper()
    mapping = next(mapper.suggest_mappings(molA, molB))

    # let's do some tests
    s = scorer(mapping)

    assert isinstance(s, float)
    np.testing.assert_allclose(actual=s, desired=0.37454, rtol=0.01)
    assert s == scorer.vals[0]
    assert scorer.i == 1

    for i in range(scorer.i, 12):
        s = scorer(mapping)
        assert s == scorer.vals[i % n_scores]

    assert scorer.i == 2


def test_build_random_dataset():
    n_compounds = 30
    compounds, mapper, scorer = build_random_dataset(n_compounds=n_compounds,
                                                     rand_seed=42)

    assert len(compounds) == n_compounds
    assert all(isinstance(c, SmallMoleculeComponent) for c in compounds)
    assert isinstance(mapper, AtomMapper)
    # assert isinstance(scorer, AtomMappingScorer)
    assert len(scorer.vals) == n_compounds


def test_build_random_mst_network():
    n_compounds = 30
    mst_network = build_random_mst_network(n_compounds=n_compounds,
                                           rand_seed=42)

    assert isinstance(mst_network, LigandNetwork)
    assert len(mst_network.nodes) == n_compounds
    assert len(mst_network.edges) == n_compounds - 1


@pytest.mark.parametrize("overlap", [0, 1, 2, 3])
@pytest.mark.parametrize("n_sub_networks", [2, 3])
def test_build_n_random_mst_network(n_sub_networks, overlap):
    n_compounds = 30
    mst_networks = build_n_random_mst_network(n_compounds=n_compounds,
                                              sub_networks=n_sub_networks,
                                              overlap=overlap, rand_seed=42)

    assert len(mst_networks) == n_sub_networks
    assert all(
        isinstance(mst_network, LigandNetwork) for mst_network in mst_networks)
    assert len(set([n for mst_network in mst_networks for n in
                    mst_network.nodes])) == n_compounds
    assert all(
        len(mst_network.edges) == len(mst_network.nodes) - 1 for mst_network in
        mst_networks)

    if overlap > 0:
        for i, net in enumerate(mst_networks):
            if i == 0:
                assert len(net.nodes) == n_compounds // n_sub_networks
            elif i == n_sub_networks - 1:
                assert len(
                    net.nodes) == n_compounds // n_sub_networks + n_compounds % n_sub_networks + overlap
            else:
                assert len(net.nodes) == n_compounds // n_sub_networks + overlap


def build_random_fully_connected_network():
    n_compounds = 30
    mst_network = build_random_fully_connected_network(n_compounds=n_compounds,
                                                       rand_seed=42)

    assert isinstance(mst_network, LigandNetwork)
    assert len(mst_network.nodes) == n_compounds
    assert len(mst_network.edges) == n_compounds * (n_compounds - 1)
