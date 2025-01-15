import random

import numpy as np

# Test Data:
nodes = range(10)
edges = [(x, y) for i, x in enumerate(nodes) for y in nodes[i + 1 :]]
weights = np.random.random_sample(len(edges))

from konnektor.network_generator_algorithms import CyclicNetworkGenerator

# Settings
# cyclesize = list(range(3,5+1))
# cyclesize = len(nodes)-1 # try to keep small!
sub_cycle_size_range = [3, 4]
node_cycle_connectivity = 2  # a node should be at least in n cycles
network_planner = CyclicNetworkGenerator(
    node_cycle_connectivity=node_cycle_connectivity, sub_cycle_size_range=sub_cycle_size_range
)
cg = network_planner.generate_network(edges=edges, weights=weights)
# %%

# Benchmark - Hyper Param test
from tqdm import tqdm

sub_cycle_size_ranges = range(3, 6)
node_cycle_connectivitys = range(1, 4)

gs = []
for sub_cycle_size_range in tqdm(sub_cycle_size_ranges, desc="cycle size"):
    gs.append(network_planner.orig_g)
    for node_cycle_connectivity in tqdm(node_cycle_connectivitys, desc="node connect", leave=False):
        print("cycleS", sub_cycle_size_range, "node con", node_cycle_connectivity)
        network_planner = CyclicNetworkGenerator(
            node_cycle_connectivity=node_cycle_connectivity,
            sub_cycle_size_range=sub_cycle_size_range,
        )
        gs.append(network_planner.generate_network(edges=edges, weights=weights))


import networkx as nx
import numpy as np
from matplotlib import pyplot as plt

cols = len(node_cycle_connectivitys) + 1
rows = len(sub_cycle_size_ranges)
fig, axes = plt.subplots(ncols=cols, nrows=rows, figsize=[16, 9])

print(gs)
gi = 0
for i, rax in enumerate(axes):
    for j, ax in enumerate(rax):
        g = gs[gi]
        print(g, ax)
        nx.draw_networkx(g, with_labels=True, ax=ax)
        if j == 0:
            ax.set_title("fully connected graph")
        else:
            ax.set_title(
                "cycles_per_node "
                + str(node_cycle_connectivitys[j - 1])
                + ", cycle_size "
                + str(sub_cycle_size_ranges[i])
            )
        gi += 1


plt.savefig("param_graphs.png")
