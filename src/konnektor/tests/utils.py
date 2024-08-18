# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest

from .data.conf import nine_mols_edges_two_sets


def test_network_merging(nine_mols_edges_two_sets):

    network1, network2 = nine_mols_edges_two_sets
