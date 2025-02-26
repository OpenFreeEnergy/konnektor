# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

from konnektor.utils.toy_data import build_random_mst_network
from konnektor.visualization import draw_ligand_network, draw_network_widget


def test_draw_ligand_network():
    """smoke test, only checking if code runs"""
    network = build_random_mst_network()
    draw_ligand_network(network)


def test_draw_network_widget():
    """smoke test, only checking if code runs"""
    network = build_random_mst_network(n_compounds=5)
    draw_network_widget(network, show_mappings=False, show_molecules=True)
