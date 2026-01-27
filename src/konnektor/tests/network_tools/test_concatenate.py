import pytest

from konnektor.network_analysis import get_is_connected
from konnektor.network_planners import MstConcatenator
from konnektor.network_tools.network_handling.concatenate import (
    append_component,
    concatenate_networks,
)
from konnektor.utils.toy_data import (
    build_n_random_mst_network,
    build_random_dataset,
    build_random_mst_network,
    emptyMapper,
    randomScorer,
)


@pytest.mark.parametrize("n_sub_networks", [2, 3, 4])
def test_concatenate_mst_networks(n_sub_networks):
    n_connecting_edges = 1
    n_compounds = 20

    networks = build_n_random_mst_network(
        n_compounds=n_compounds, sub_networks=n_sub_networks, overlap=0, rand_seed=42
    )
    concatenator = MstConcatenator(
        emptyMapper(),
        randomScorer(n_scores=n_compounds**2),
        n_connecting_edges=n_connecting_edges,
    )

    new_network = concatenate_networks(networks, concatenator)

    assert len(networks) == n_sub_networks
    assert len(new_network.nodes) == n_compounds
    # network edges + the network connecting edges
    assert len(new_network.edges) == sum(
        [len(n.edges) for n in networks]
    ) + n_connecting_edges * sum([i for i in range(1, n_sub_networks)])
    assert get_is_connected(new_network)


def test_append_node():
    n_connecting_edges = 2
    n_compounds = 20

    network = build_random_mst_network(n_compounds=n_compounds - 1, rand_seed=42)
    compounds, _, _ = build_random_dataset(n_compounds=1, rand_seed=42)
    concatenator = MstConcatenator(
        emptyMapper(),
        randomScorer(n_scores=n_compounds**2),
        n_connecting_edges=n_connecting_edges,
    )

    new_network = append_component(network, compounds[0], concatenator=concatenator)

    assert len(new_network.nodes) == n_compounds
    # network edges + the network connecting edges
    assert len(new_network.edges) == len(network.edges) + n_connecting_edges
    assert get_is_connected(new_network)
