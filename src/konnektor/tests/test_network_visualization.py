# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

import konnektor
from konnektor.tests.network_planners.conf import GenAtomMapper
from konnektor.utils.toy_data import build_random_mst_network
from konnektor.visualization import draw_ligand_network, draw_network_widget


@pytest.fixture
def toluene_star_network(toluene_vs_others):
    toluene, others = toluene_vs_others
    mapper = GenAtomMapper()

    planner = konnektor.network_planners.RadialNetworkGenerator(mappers=mapper, scorer=None)
    return planner.generate_ligand_network(components=others, central_component=toluene)


def test_draw_ligand_network():
    """smoke test, only checking if code runs"""
    network = build_random_mst_network()
    draw_ligand_network(network)


def test_draw_network_widget():
    """smoke test, only checking if code runs"""
    network = build_random_mst_network(n_compounds=5)
    draw_network_widget(network, show_mappings=False, show_molecules=True)


def test_draw_star_ligand_network(toluene_star_network):
    """smoke test, only checking if code runs"""
    draw_ligand_network(toluene_star_network)
