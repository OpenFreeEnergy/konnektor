# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import numpy as np

from konnektor.network_analysis import (
    get_component_number_cycles,
    get_is_connected,
    get_network_score,
)
from konnektor.network_planners import CyclicNetworkGenerator
from konnektor.utils.toy_data import build_random_dataset


def test_cyclic_network_planner():
    n_compounds = 8
    ncycles = 2
    components, emptyMapper, genScorer = build_random_dataset(n_compounds=n_compounds, rand_seed=42)

    planner = CyclicNetworkGenerator(
        mappers=emptyMapper,
        scorer=genScorer,
        cycle_sizes=3,
        node_present_in_cycles=ncycles,
    )
    network = planner.generate_ligand_network(components)

    assert len(network.nodes) == n_compounds
    edge_count = n_compounds * 3
    assert len(network.edges) <= edge_count
    assert len(network.edges) > n_compounds
    assert get_is_connected(network)
    nnode_cycles = get_component_number_cycles(network)
    assert all(v >= ncycles for k, v in nnode_cycles.items())

    np.testing.assert_allclose(get_network_score(network), 10.347529, rtol=0.01)
