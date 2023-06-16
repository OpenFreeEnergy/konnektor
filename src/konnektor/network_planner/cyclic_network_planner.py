import itertools
import logging
import numpy as np
from datetime import datetime

from typing import Iterable, List, Tuple

import networkx as nx
from networkx import Graph

from gufe import SmallMoleculeComponent, AtomMapper

from openfe.setup.ligand_network import LigandNetwork  # only temproary
from ._abstract_network_planner import _AbstractNetworkPlanner, Network

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


class cyclicNetwork(_AbstractNetworkPlanner):

    def __init__(self, node_cycle_connectivity: int = 2,
                 sub_cycle_size_range: Iterable[int] = 3):
        self.node_cycle_connectivity = node_cycle_connectivity

        if (isinstance(sub_cycle_size_range, int)):
            self.max_sub_cycle_size = sub_cycle_size_range
            self.sub_cycle_size_range = [sub_cycle_size_range]
        elif (isinstance(sub_cycle_size_range, Iterable) and all([isinstance(i, int) for i in sub_cycle_size_range])):
            self.max_sub_cycle_size = max(sub_cycle_size_range)
            self.sub_cycle_size_range = sub_cycle_size_range
        else:
            raise ValueError("Cycle needs to be int or list of ints")

        if (self.max_sub_cycle_size < 3):
            raise ValueError("a cycle has a minimal size of 3.")

        self.orig_g = None
        self.cycle_metric = self._score_summation  # how to evaluate the cost of a cycle.

    def generate_network(self, edges, weights) -> Network:
        pass

    def generate_cyclic_graph(self, edges: List[Tuple[int, int]], weights: List[float]) -> Graph:
        log.info("Building Cyclic Graph - START")
        start_time_total = datetime.now()

        self.orig_g = self._translate_input(edges, weights)

        # Generate paths of given size
        start_time_path_generation = datetime.now()
        log.info("Cycle Selection init - " + str(start_time_path_generation))

        self._cyclic_paths_scored = self._generate_cyclic_paths(graph=self.orig_g)

        end_time_path_generation = datetime.now()
        duration_path_generation = end_time_path_generation - start_time_path_generation
        log.info("Duration: " + str(duration_path_generation))
        log.info("Cycle Selection complete")

        # Select Cycles, to give each node a given connectivity:
        start_time_cycle_selection = datetime.now()
        log.info("Cycle Selection init - " + str(start_time_cycle_selection))

        self._selected_cycles = self._greedy_select_cylces_by_node_connectivity(
            cyclic_paths_scored=self._cyclic_paths_scored, graph=self.orig_g)

        end_time_cycle_selection = datetime.now()
        duration_cycle_selection = end_time_cycle_selection - start_time_cycle_selection
        log.info("Duration: " + str(duration_cycle_selection))
        log.info("Cycle Selection complete")

        cyclic_graph = self._reduce_graph_to_cycles(selected_cycles=self._selected_cycles, graph=self.orig_g)

        log.info("Building Cyclic Graph - DONE")
        end_time_total = datetime.now()
        duration_total = end_time_total - start_time_total

        log.info("Timings:")
        log.info("\t Cycle generation duration: " + str(duration_path_generation))
        log.info("\t Cycle selection duration: " + str(duration_cycle_selection))
        log.info("\t total duration: " + str(duration_total))

        return cyclic_graph

    """
    Cycle building
    """

    def _translate_input(self, edges: List[Tuple[int, int]], weights: List[float]) -> Graph:
        log.info("Translate input to Networkx Graph")

        # build Edges:
        w_edges = []
        nodes = []
        for e, w in zip(edges, weights):
            w_edges.append((e[0], e[1], w))
            nodes.extend(e)
        nodes = list(set(nodes))

        # Build graph
        g = nx.Graph()
        [g.add_node(n) for n in nodes]
        g.add_weighted_edges_from(ebunch_to_add=w_edges)

        return g

    def _reduce_graph_to_cycles(self, selected_cycles, graph: Graph) -> Graph:
        # Build new Edges for New Graph
        log.info("Build Graph with new edge selection")
        new_edges = []
        for c in selected_cycles:
            for i in range(len(c)):
                n1 = c[i]
                n2 = c[(i + 1) % len(c)]
                old_e = graph.get_edge_data(n1, n2)
                s = old_e['weight']
                new_edges.append((n1, n2, s))
        new_edges = list(set(new_edges))

        # Build graph
        g = nx.Graph()
        [g.add_node(n) for n in graph.nodes]
        g.add_weighted_edges_from(ebunch_to_add=new_edges)

        return g

    def _generate_cyclic_paths(self, graph: Graph) -> Tuple[List[int], float]:
        n_nodes = len(graph.nodes)
        if (self.max_sub_cycle_size > n_nodes):
            raise ValueError("a cycle can not be larger, than the number of nodes.")

        # Build Paths
        log.info("Generate all paths up to size " + str(self.max_sub_cycle_size))
        raw_cycles = [c for c in nx.simple_cycles(graph, length_bound=self.max_sub_cycle_size)]
        log.info("Found  " + str(len(raw_cycles)) + " cycles")

        # Filter cycle sizes:
        log.info("Filter too small cycles")
        selected_cycle_by_size = set([tuple(x) for x in raw_cycles if (len(x) in self.sub_cycle_size_range)])

        # Get Cycle Weigths -> Dist matrix
        log.info("Aggregate edge weight to cycle weigth")
        cyclic_paths_scored = []
        for c in selected_cycle_by_size:
            d = self.cycle_metric(list(c), graph)  # cycle cost
            cyclic_paths_scored.append((c, np.round(d, 2)))
        cyclic_paths_scored = np.array(cyclic_paths_scored, dtype=object)

        return cyclic_paths_scored

    def _greedy_select_cylces_by_node_connectivity(self, cyclic_paths_scored: Tuple[List[int], float], graph: Graph) -> \
    List[List[int]]:
        # Build priorityQueue for MST selection
        log.info("Build Cycle Priority Queue")
        cycle_priority_queue = list(map(lambda x: x[0], sorted(cyclic_paths_scored, key=lambda x: x[1])))

        # Greedy opt algorithm: optimize connectivity criteria + select cycles Kruskal-MST like
        log.info("Optimize such that each node appears in " + str(self.node_cycle_connectivity) + " Cycles")
        connectivity_dict = {k: 0 for k in graph}
        selected_cycles = []

        # opt criteria
        termination_criteria = lambda d: all([v >= self.node_cycle_connectivity for k, v in d.items()])
        while not termination_criteria(connectivity_dict):

            # MST like iteration
            selected_cycles = []
            for c in cycle_priority_queue:
                if termination_criteria(connectivity_dict):  # efficiency break
                    break
                elif (any([connectivity_dict[n] < self.node_cycle_connectivity for n in c])):
                    selected_cycles.append(c)
                    for n in c:
                        connectivity_dict[n] += 1
                else:
                    continue

        log.info("Number  of required Cycles: " + str(len(selected_cycles)))
        log.info("Node Cycle appearance: " + str(connectivity_dict))

        return selected_cycles

    """
        Cycle metrics
    """

    def _score_summation(self, c: List[int], g: Graph) -> float:
        d = 0
        for i in range(len(c)):  # Calculate the score of the cycle
            n1 = c[i]
            n2 = c[(i + 1) % len(c)]
            e = g.get_edge_data(n1, n2)
            d += e['weight']

        return d
