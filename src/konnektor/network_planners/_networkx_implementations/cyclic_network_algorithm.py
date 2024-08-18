# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import logging
from datetime import datetime
from typing import Iterable, List, Tuple

import networkx as nx
import numpy as np
from networkx import Graph

from ._abstract_network_algorithm import _AbstractNetworkAlgorithm

log = logging.getLogger(__name__)
log.setLevel(logging.WARNING)


class CyclicNetworkAlgorithm(_AbstractNetworkAlgorithm):

    def __init__(
        self, node_cycle_connectivity: int = 2, sub_cycle_size_range: Iterable[int] = 3
    ):
        self.node_cycle_connectivity = node_cycle_connectivity

        if isinstance(sub_cycle_size_range, int):
            self.max_sub_cycle_size = sub_cycle_size_range
            self.sub_cycle_size_range = [sub_cycle_size_range]
        elif isinstance(sub_cycle_size_range, Iterable) and all(
            [isinstance(i, int) for i in sub_cycle_size_range]
        ):
            self.max_sub_cycle_size = max(sub_cycle_size_range)
            self.sub_cycle_size_range = sub_cycle_size_range
        else:
            raise ValueError("Cycle needs to be int or list of ints")

        if self.max_sub_cycle_size < 3:
            raise ValueError("a cycle has a minimal size of 3.")

        self.orig_g = None
        self.cycle_metric = (
            self._score_summation
        )  # how to evaluate the cost of a cycle.

    def generate_network(
        self, edges: list[tuple[int, int]], weights: list[float]
    ) -> Graph:
        log.info("Building Cyclic Graph - START")
        start_time_total = datetime.now()

        self.orig_g = self._translate_input(edges, weights)

        # Generate paths of given size
        start_time_path_generation = datetime.now()
        log.info("Cycle Selection init - " + str(start_time_path_generation))

        self._cyclic_paths_scored = self._generate_cyclic_paths(graph=self.orig_g)

        end_time_path_generation = datetime.now()
        duration_path_generation = end_time_path_generation - start_time_path_generation
        log.info("\tDuration: " + str(duration_path_generation))
        log.info("Cycle Selection complete\n")

        # Select Cycles, to give each node a given connectivity:
        start_time_cycle_selection = datetime.now()
        log.info("Cycle Selection init - " + str(start_time_cycle_selection))

        self._selected_cycles = self._greedy_select_cylces_by_node_connectivity(
            cyclic_paths_scored=self._cyclic_paths_scored, graph=self.orig_g
        )

        end_time_cycle_selection = datetime.now()
        duration_cycle_selection = end_time_cycle_selection - start_time_cycle_selection
        log.info("\tDuration: " + str(duration_cycle_selection))
        log.info("Cycle Selection complete\n")

        self.cyclic_graph = self._reduce_graph_to_cycles(
            selected_cycles=self._selected_cycles, graph=self.orig_g
        )

        log.info("Building Cyclic Graph - DONE")
        end_time_total = datetime.now()
        duration_total = end_time_total - start_time_total

        log.info("TIMINGS:")
        log.info("\t Cycle generation duration: " + str(duration_path_generation))
        log.info("\t Cycle selection duration: " + str(duration_cycle_selection))
        log.info("\t total duration: " + str(duration_total))

        return self.cyclic_graph

    """
    Cycle building
    """

    def _translate_input(
        self, edges: List[Tuple[int, int]], weights: List[float]
    ) -> Graph:
        log.info("\tTranslate input to Networkx Graph")

        # build Edges:
        w_edges = []
        nodes = []
        # The initial "weights" are Scores, which need to be translated to weights.
        weights = list(map(lambda x: 1 - x, weights))
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
        log.info("\tBuild Graph with new edge selection")
        new_edges = []
        for c in selected_cycles:
            for i in range(len(c)):
                n1 = c[i]
                n2 = c[(i + 1) % len(c)]
                old_e = graph.get_edge_data(n1, n2)
                s = old_e["weight"]
                new_edges.append((n1, n2, s))
        new_edges = list(set(new_edges))

        # Build graph
        g = nx.Graph()
        [g.add_node(n) for n in graph.nodes]
        g.add_weighted_edges_from(ebunch_to_add=new_edges)

        return g

    def _generate_cyclic_paths(self, graph: Graph) -> Tuple[List[int], float]:
        n_nodes = len(graph.nodes)
        if self.max_sub_cycle_size > n_nodes:
            raise ValueError("a cycle can not be larger, than the number of nodes.")

        # Build Paths
        log.info("\tGenerate all paths up to size " + str(self.max_sub_cycle_size))
        raw_cycles = [
            c for c in nx.simple_cycles(graph, length_bound=self.max_sub_cycle_size)
        ]
        log.info("\tFound  " + str(len(raw_cycles)) + " cycles")

        # Filter cycle sizes:
        log.info("\tFilter too small cycles")
        selected_cycle_by_size = set(
            [tuple(x) for x in raw_cycles if (len(x) in self.sub_cycle_size_range)]
        )

        # Get Cycle Weigths -> Dist matrix
        log.info("\tAggregate edge weight to cycle weigth")
        cyclic_paths_scored = []
        for c in selected_cycle_by_size:
            d = self.cycle_metric(list(c), graph)  # cycle cost
            cyclic_paths_scored.append((c, np.round(d, 2)))
        cyclic_paths_scored = np.array(cyclic_paths_scored, dtype=object)

        return cyclic_paths_scored

    def _greedy_select_cylces_by_node_connectivity(
        self, cyclic_paths_scored: Tuple[List[int], float], graph: Graph
    ) -> List[List[int]]:
        # Build priorityQueue for MST selection
        log.info("\tBuild Cycle Priority Queue")
        cycle_priority_queue = list(
            map(lambda x: x[0], sorted(cyclic_paths_scored, key=lambda x: x[1]))
        )

        # Greedy opt algorithm: optimize connectivity criteria + select cycles Kruskal-MST like
        log.info(
            "\tOptimize such that each node appears in "
            + str(self.node_cycle_connectivity)
            + " Cycles"
        )
        cycle_connectivity_dict = {k: 0 for k in graph}

        # opt criteria
        i = 0
        j = 0
        termination_criteria = lambda d: all(
            [v >= self.node_cycle_connectivity for k, v in d.items()]
        )
        selected_cycles = []
        # print("cycles", len(cycle_priority_queue))
        while not termination_criteria(cycle_connectivity_dict):
            # print("ITER", j, i)

            # MST like iteration
            c = cycle_priority_queue[i]
            if (
                any(
                    [
                        cycle_connectivity_dict[n] < self.node_cycle_connectivity
                        for n in c
                    ]
                )
                and c not in selected_cycles
            ):
                # print("select cycle:", c)
                selected_cycles.append(c)
                for n in c:
                    cycle_connectivity_dict[n] += 1

                # Be efficient, remove already counted cycles
                cyclic_paths_scored = list(
                    filter(lambda x: x[0] != c, cyclic_paths_scored)
                )

                # update cycle score - punish by note appearance:
                resort_f = lambda x: x[1] * max(
                    [cycle_connectivity_dict[i] for i in x[0]]
                )
                cycle_priority_queue = list(
                    map(lambda x: x[0], sorted(cyclic_paths_scored, key=resort_f))
                )
                i = 0
            else:
                i = i + 1 % len(cycle_priority_queue)

            # no inf loop!
            if j > 1000:
                # print("did not converge!")
                break
            else:
                j += 1

        # print(cycle_connectivity_dict)
        self.cycle_connectivity_dict = cycle_connectivity_dict
        log.info("\tNumber  of required Cycles: " + str(len(selected_cycles)))
        log.info("\tNode Cycle appearance: " + str(cycle_connectivity_dict))

        return selected_cycles

    def generate_network_double_greedy(
        self, edges: List[Tuple[int, int]], weights: List[float], edge_limitor=200
    ) -> Graph:
        log.info("Building Cyclic Graph - START")
        start_time_total = datetime.now()

        ew = dict(zip(map(lambda x: tuple(sorted(list(x))), edges), weights))
        nodes = [n for e in edges for n in e]

        start_time_cycle_generation = datetime.now()
        log.info("Priority Queue Gen init - " + str(start_time_cycle_generation))

        # build
        self.orig_g = self._translate_input(edges=edges, weights=weights)
        from collections import defaultdict

        min_node = defaultdict(list)
        for (e1, e2), v in sorted(ew.items(), key=lambda x: x[1]):
            if len(min_node) == len(nodes) and not all(
                [len(min_node[n]) < edge_limitor for n in min_node]
            ):
                break
            if len(min_node[e1]) < edge_limitor:
                min_node[e1].append((e1, e2))
            if len(min_node[e2]) < edge_limitor:
                min_node[e2].append((e1, e2))

        end_time_cycle_generation = datetime.now()
        duration_cycle_generation = (
            end_time_cycle_generation - start_time_cycle_generation
        )
        log.info("\tDuration: " + str(duration_cycle_generation))
        log.info("Priority Queue Gen complete\n")

        start_time_cycle_selection = datetime.now()
        log.info("Cycle Selection init - " + str(start_time_cycle_selection))

        # build_cycles for each node:
        all_cycles_per_node = {}
        all_edge_collection = []
        all_node_collection = []
        for target_node, prioqueue in min_node.items():
            not_target_node = lambda e: not (e[0] == target_node or e[1] == target_node)

            node_cycles = []
            for i in range(self.node_cycle_connectivity):
                start_edge = prioqueue[i]
                cycle_edges = [start_edge]  # start edge
                cycle_nodes = set(start_edge)
                all_edge_collection.append(start_edge)
                all_node_collection.extend(list(start_edge))

                adding_node = lambda e: len(cycle_nodes.intersection(set(e))) == 1
                use_edge = (
                    lambda e: not e in cycle_edges
                    and adding_node(e)
                    and not_target_node(e)
                )

                successor_node = (
                    prioqueue[i][1]
                    if (target_node == prioqueue[i][0])
                    else prioqueue[i][0]
                )
                for j in range(self.max_sub_cycle_size - 2):
                    tmp_prio_queue = min_node[successor_node]

                    min_es = [e for e in tmp_prio_queue if (use_edge(e))]
                    if len(min_es) == 0:
                        break
                    else:
                        min_e = min_es[0]

                    cycle_edges.append(min_e)

                    cycle_nodes = cycle_nodes.union(min_e)
                    all_edge_collection.append(tuple(min_e))
                    successor_node = (
                        min_e[1] if (successor_node == min_e[0]) else min_e[0]
                    )

                if len(cycle_edges) == self.max_sub_cycle_size - 1:
                    cycle_edges.append(tuple(sorted([successor_node, target_node])))
                    node_cycles.append(cycle_edges)
            all_cycles_per_node[target_node] = node_cycles

        self.all_cycles = set(
            [tuple(c) for cs in all_cycles_per_node.values() for c in cs]
        )

        all_edges = set([e for c in self.all_cycles for e in c])
        end_time_cycle_selection = datetime.now()
        duration_cycle_selection = end_time_cycle_selection - start_time_cycle_selection
        log.info("\tDuration: " + str(duration_cycle_selection))
        log.info("Cycle Selection complete\n")

        log.info("Building Cyclic Graph - DONE")
        end_time_total = datetime.now()
        duration_total = end_time_total - start_time_total

        log.info("TIMINGS:")
        log.info("\t Cycle generation duration: " + str(duration_cycle_generation))
        log.info("\t Cycle selection duration: " + str(duration_cycle_selection))
        log.info("\t total duration: " + str(duration_total))

        self.cyclic_graph = nx.Graph()
        [self.cyclic_graph.add_node(n) for n in nodes]
        self.cyclic_graph.add_weighted_edges_from(
            ebunch_to_add=[(e[0], e[1], ew[e]) for e in all_edges]
        )
        return self.cyclic_graph

    """
        Cycle metrics
    """

    def _score_summation(self, c: List[int], g: Graph) -> float:
        d = 0
        for i in range(len(c)):  # Calculate the score of the cycle
            n1 = c[i]
            n2 = c[(i + 1) % len(c)]
            e = g.get_edge_data(n1, n2)
            d += e["weight"]

        return d
