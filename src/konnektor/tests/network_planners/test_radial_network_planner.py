# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

import konnektor
from konnektor import network_planners

from .conf import (atom_mapping_basic_test_files, toluene_vs_others,
                   mol_from_smiles, genScorer,
                   GenAtomMapper, BadMapper, ErrorMapper)

from gufe import SmallMoleculeComponent


@pytest.mark.parametrize('as_list', [False])
def test_radial_network(atom_mapping_basic_test_files, toluene_vs_others,
                        as_list):
    toluene, others = toluene_vs_others
    central_ligand_name = 'toluene'
    mapper = GenAtomMapper()

    if as_list:
        mapper = [mapper]

    planner = konnektor.network_planners.RadialLigandNetworkPlanner(mapper=mapper, scorer=None)
    network = planner.generate_ligand_network(ligands=others, central_ligand=toluene)

    # couple sanity checks
    assert len(network.nodes) == len(atom_mapping_basic_test_files)
    assert len(network.edges) == len(others)
    # check that all ligands are present, i.e. we included everyone
    ligands_in_network = {mol.name for mol in network.nodes}
    assert ligands_in_network == set(atom_mapping_basic_test_files.keys())
    # check that every edge has the central ligand within
    assert all((central_ligand_name in {mapping.componentA.name, mapping.componentB.name})
               for mapping in network.edges)


def test_radial_network_with_scorer(toluene_vs_others):
    toluene, others = toluene_vs_others

    mapper = GenAtomMapper()
    scorer = genScorer
    planner = konnektor.network_planners.RadialLigandNetworkPlanner(mapper=mapper, scorer=scorer)
    network = planner.generate_ligand_network(ligands=others, central_ligand=toluene)

    assert len(network.edges) == len(others)

    for edge in network.edges:
        assert len(edge.componentA_to_componentB) > 1  # we didn't take the bad mapper
        assert 'score' in edge.annotations
        assert edge.annotations['score'] == 1.0 / len(edge.componentA_to_componentB)


def test_radial_network_multiple_mappers_no_scorer(toluene_vs_others):
    toluene, others = toluene_vs_others
    # in this one, we should always take the bad mapper
    mapper = BadMapper()
    planner = konnektor.network_planners.RadialLigandNetworkPlanner(mapper=mapper, scorer=None)
    network = planner.generate_ligand_network(ligands=others, central_ligand=toluene)

    assert len(network.edges) == len(others)

    for edge in network.edges:
        assert edge.componentA_to_componentB == {0: 0}


def test_radial_network_failure(atom_mapping_basic_test_files):
    nigel = SmallMoleculeComponent(mol_from_smiles('N'))

    mapper = ErrorMapper()
    planner = konnektor.network_planners.RadialLigandNetworkPlanner(mapper=mapper, scorer=None)

    with pytest.raises(ValueError, match='No mapping found for'):
        network = planner.generate_ligand_network(ligands=[nigel], central_ligand=atom_mapping_basic_test_files['toluene'])

