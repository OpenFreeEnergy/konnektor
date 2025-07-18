import datetime as dt

import numpy as np

from konnektor.network_generator_algorithms import CyclicNetworkGenerator

sub_cycle_size_range = 3
node_cycle_connectivity = 2  # a node should be at least in n cycles

# Test Data: $10, 20, 30, 40, 50, 60,70,80, 90, 100, 300, 500, 700, 1000, 5000, 6000,
steps = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 300]  #
cgs = []
for n in steps:  # range(10, 81, 10):
    print("Number of nodes:", n)
    nodes = range(n)
    edges = [(x, y) for i, x in enumerate(nodes) for y in nodes[i + 1 :]]
    weights = np.random.random_sample(len(edges))

    network_planner = CyclicNetworkGenerator(
        node_cycle_connectivity=node_cycle_connectivity, sub_cycle_size_range=sub_cycle_size_range
    )
    start = dt.datetime.now()
    cg = network_planner.generate_network(edges=edges, weights=weights)
    cgs.append([network_planner.orig_g, cg])
    end = dt.datetime.now()
    duration = end - start
    print("  Duration: ", duration)
    print()
