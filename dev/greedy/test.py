import random
import numpy as np

#Test Data:
nodes = (range(10))
edges = [(x,y) for i, x in enumerate(nodes) for y in nodes[i+1:]]
weights = np.random.random_sample(len(edges))

from konnektor.network_generator_algorithms import CyclicNetworkGenerator

# Settings
#cyclesize = list(range(3,5+1))
#cyclesize = len(nodes)-1 # try to keep small!
sub_cycle_size_range = [3,4,5]

node_cycle_connectivity = 3 #a node should be at least in n cycles
network_planner = CyclicNetworkGenerator(node_cycle_connectivity=node_cycle_connectivity, sub_cycle_size_range=sub_cycle_size_range)

cg = network_planner.generate_network(edges=edges, weights=weights)


#Plooot
import networkx as nx
from matplotlib import pyplot as plt

fig, axes = plt.subplots(ncols=2, nrows=1, figsize=[16,9])

nx.draw_networkx(network_planner.orig_g, with_labels=True, ax=axes[0])
nx.draw_networkx(cg, with_labels=True, ax=axes[1])

axes[0].set_title("fully connected graph")
axes[1].set_title("cycles_per_node "+str(network_planner.node_cycle_connectivity)+", cycle_size "+str(network_planner.sub_cycle_size_range))

plt.savefig("cycGraph.png")