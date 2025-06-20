# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
from gufe import SmallMoleculeComponent

import konnektor
from konnektor.tests.network_planners.conf import (
    BadMapper,
    BadMultiMapper,
    ErrorMapper,
    GenAtomMapper,
    SuperBadMapper,
    genScorer,
    mol_from_smiles,
)


@pytest.mark.parametrize("as_list", [False, True])
def test_star_network(atom_mapping_basic_test_files, toluene_vs_others, as_list):
    toluene, others = toluene_vs_others
    central_ligand_name = "toluene"
    mapper = GenAtomMapper()

    if as_list:
        mapper = [mapper]

    planner = konnektor.network_planners.RadialNetworkGenerator(mappers=mapper, scorer=None)
    network = planner.generate_ligand_network(components=others, central_component=toluene)

    # couple sanity checks
    assert len(network.nodes) == len(atom_mapping_basic_test_files)
    assert len(network.edges) == len(others)
    # check that all ligands are present, i.e. we included everyone
    ligands_in_network = {mol.name for mol in network.nodes}
    assert ligands_in_network == set(atom_mapping_basic_test_files.keys())
    # check that every edge has the central ligand within
    assert all(
        (central_ligand_name in {mapping.componentA.name, mapping.componentB.name})
        for mapping in network.edges
    )


def test_star_network_central_in_components(atom_mapping_basic_test_files, toluene_vs_others):
    """issue #544, if including the central ligand in "ligands", warn and make sure there isn't self-edge"""
    toluene, others = toluene_vs_others
    all_ligands = others + [toluene]
    central_ligand_name = "toluene"
    mapper = GenAtomMapper()

    planner = konnektor.network_planners.RadialNetworkGenerator(mappers=mapper, scorer=None)
    with pytest.warns(UserWarning, match="The central component 'toluene' is present in the list"):
        network = planner.generate_ligand_network(components=all_ligands, central_component=toluene)

    assert ("toluene", "toluene") not in {
        (e.componentA.name, e.componentB.name) for e in network.edges
    }
    # couple sanity checks
    assert len(network.nodes) == len(atom_mapping_basic_test_files)
    assert len(network.edges) == len(others)
    # check that all ligands are present, i.e. we included everyone
    ligands_in_network = {mol.name for mol in network.nodes}
    assert ligands_in_network == set(atom_mapping_basic_test_files.keys())
    # check that every edge has the central ligand within
    assert all(
        (central_ligand_name in {mapping.componentA.name, mapping.componentB.name})
        for mapping in network.edges
    )


def test_star_network_with_scorer(toluene_vs_others):
    toluene, others = toluene_vs_others

    mappers = [BadMapper(), GenAtomMapper()]
    scorer = genScorer
    planner = konnektor.network_planners.RadialNetworkGenerator(mappers=mappers, scorer=scorer)
    network = planner.generate_ligand_network(components=others, central_component=toluene)

    assert len(network.edges) == len(others)

    for edge in network.edges:
        assert len(edge.componentA_to_componentB) > 1  # we didn't take the bad mapper
        assert "score" in edge.annotations
        assert edge.annotations["score"] == (1 - 1.0 / len(edge.componentA_to_componentB))


def test_star_network_multiple_mappers_no_scorer(toluene_vs_others):
    toluene, others = toluene_vs_others
    # in this one, we should always take the bad mapper
    mappers = [BadMultiMapper(), SuperBadMapper(), BadMapper()]
    planner = konnektor.network_planners.RadialNetworkGenerator(mappers=mappers, scorer=None)

    with pytest.warns(match="Only the first valid mapper will be used: <BadMulti"):
        network = planner.generate_ligand_network(components=others, central_component=toluene)

    assert len(network.edges) == len(others)
    for edge in network.edges:
        assert edge.componentA_to_componentB == {0: 2}


def test_star_network_multiple_mappers_no_scorer_takes_first_valid(toluene_vs_others):
    toluene, others = toluene_vs_others
    # in this one, we should always take the bad mapper
    mappers = [ErrorMapper(), BadMultiMapper()]
    planner = konnektor.network_planners.RadialNetworkGenerator(mappers=mappers, scorer=None)

    with pytest.warns(match="Only the first valid mapper will be used: <BadMulti"):
        network = planner.generate_ligand_network(components=others, central_component=toluene)

    assert len(network.edges) == len(others)
    for edge in network.edges:
        assert edge.componentA_to_componentB == {0: 2}


def test_star_network_failure(atom_mapping_basic_test_files):
    nigel = SmallMoleculeComponent(mol_from_smiles("N"))

    mapper = ErrorMapper()
    planner = konnektor.network_planners.RadialNetworkGenerator(mappers=mapper, scorer=None)

    with pytest.raises(RuntimeError, match="Could not generate any mapping"):
        planner.generate_ligand_network(
            components=[nigel],
            central_component=atom_mapping_basic_test_files["toluene"],
        )
