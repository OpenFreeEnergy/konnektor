# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import pytest
from gufe import LigandNetwork

from konnektor.network_planners.concatenators import (
    MaxConcatenator,
    RedundantMstConcatenator,
)
from konnektor.utils.toy_data import (
    EmptyMapper,
    RandomScorer,
    build_n_random_mst_network,
    build_random_dataset,
)


def _bridges(concatenated, subnetworks):
    """Get the connecting edges (those not present in input subnetworks)."""
    original = set().union(*(n.edges for n in subnetworks))
    return [e for e in concatenated.edges if e not in original]


def _pairs(mappings):
    return {frozenset((m.componentA, m.componentB)) for m in mappings}


@pytest.fixture
def three_subnetworks():
    n = 30
    networks = build_n_random_mst_network(
        n_compounds=n,
        sub_networks=3,
        overlap=0,
        rand_seed=42,
    )
    return n, networks


@pytest.fixture
def two_singleton_subnetworks():
    (network,) = build_n_random_mst_network(
        n_compounds=2,
        sub_networks=1,
        overlap=0,
        rand_seed=1,
    )
    x, y = list(network.nodes)

    return [
        LigandNetwork(nodes=[x], edges=[]),
        LigandNetwork(nodes=[y], edges=[]),
    ]


@pytest.mark.parametrize("n_sub_networks", [2, 3, 4])
def test_redundant_overlays_n_spanning_trees(n_sub_networks):
    """n_redundancy overlaid trees add n_redundancy*(k-1) edges."""
    n = 30
    n_redundancy = 2
    networks = build_n_random_mst_network(
        n_compounds=n, sub_networks=n_sub_networks, overlap=0, rand_seed=42
    )
    concatenator = RedundantMstConcatenator(
        EmptyMapper(), RandomScorer(n=n), n_redundancy=n_redundancy
    )
    result = concatenator.concatenate_networks(networks)

    assert result.is_connected()
    assert len(result.nodes) == n
    n_new = len(_bridges(result, networks))
    assert n_new == n_redundancy * (n_sub_networks - 1)


def test_redundant_does_not_reuse_bridge_pairs(three_subnetworks):
    """Each pass excludes earlier edges, so the bridges don't repeat."""
    n, networks = three_subnetworks
    result = RedundantMstConcatenator(
        EmptyMapper(), RandomScorer(n=n), n_redundancy=2
    ).concatenate_networks(networks)

    bridges = _bridges(result, networks)
    # no bridge pair is used twice across the overlaid trees
    assert len(_pairs(bridges)) == len(bridges)


def test_redundant_raises_if_first_tree_cannot_be_built(two_singleton_subnetworks):
    """Failure to build the first spanning tree raises."""
    concatenator = RedundantMstConcatenator(
        EmptyMapper(),
        RandomScorer(n=2),
        n_redundancy=2,
    )
    network_a, network_b = two_singleton_subnetworks
    mappings = concatenator._score_pair_edges(network_a, network_b, exclude=set())
    assert len(mappings) == 1

    with pytest.raises(RuntimeError, match="Could not build"):
        concatenator.concatenate_networks(
            two_singleton_subnetworks,
            exclude_edges=[mappings[0]],
        )


def test_redundant_warns_if_later_tree_cannot_be_built(
    two_singleton_subnetworks,
):
    """Failure to build a later spanning tree warns and keeps earlier trees."""
    concatenator = RedundantMstConcatenator(
        EmptyMapper(),
        RandomScorer(n=2),
        n_redundancy=2,
    )

    with pytest.warns(UserWarning, match="Could only build 1"):
        result = concatenator.concatenate_networks(two_singleton_subnetworks)

    assert result.is_connected()
    assert len(_bridges(result, two_singleton_subnetworks)) == 1


def test_redundant_exclude_edges_never_used(three_subnetworks):
    """aexclude_edges are excluded across every pass."""
    n, networks = three_subnetworks
    concatenator = RedundantMstConcatenator(EmptyMapper(), RandomScorer(n=n), n_redundancy=2)

    base = concatenator.concatenate_networks(networks)
    exclude = [_bridges(base, networks)[0]]

    result = concatenator.concatenate_networks(networks, exclude_edges=exclude)
    assert _pairs(exclude).isdisjoint(_pairs(_bridges(result, networks)))
    assert result.is_connected()


def test_redundant_partial_forest_not_counted_as_a_tree():
    """A pass that can only form a partial forest is discarded, not counted."""
    n = 4
    components, _, _ = build_random_dataset(n_compounds=n, rand_seed=1)
    c0, c1, c2, c3 = components
    frags = [LigandNetwork(nodes=[c], edges=[]) for c in components]

    all_cross = _bridges(
        MaxConcatenator(EmptyMapper(), RandomScorer(n=n)).concatenate_networks(frags), frags
    )
    assert len(all_cross) == 6  # complete graph over 4 singletons

    # Keep exactly four edges. The first spanning tree consumes three, so only
    # one edge remains and a second complete spanning tree is impossible.
    keep = {
        frozenset((c0, c1)),
        frozenset((c1, c2)),
        frozenset((c2, c3)),
        frozenset((c1, c3)),
    }
    exclude = [e for e in all_cross if frozenset((e.componentA, e.componentB)) not in keep]
    assert len(exclude) == 2  # we kept 4 edges

    concatenator = RedundantMstConcatenator(EmptyMapper(), RandomScorer(n=n), n_redundancy=2)
    with pytest.warns(UserWarning, match="Could only build"):
        result = concatenator.concatenate_networks(frags, exclude_edges=exclude)

    # only the first full tree's bridges remain; the partial 2nd pass is discarded
    assert result.is_connected()
    assert len(_bridges(result, frags)) == n - 1


@pytest.mark.parametrize("n_redundancy", [0, -1])
def test_redundancy_must_be_positive(n_redundancy):
    with pytest.raises(ValueError, match="at least 1"):
        RedundantMstConcatenator(
            EmptyMapper(),
            RandomScorer(n=10),
            n_redundancy=n_redundancy,
        )
