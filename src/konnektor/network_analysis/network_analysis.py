import random
import numpy as np
import networkx as nx
from gufe import LigandNetwork

from .. import network_tools as tools


def get_is_connected(cg:LigandNetwork)-> bool:
    return nx.is_connected(cg.graph.to_undirected())


def get_graph_score(cg:LigandNetwork)->float:
    score = sum([e.annotations["score"] for e in cg.edges])
    return score


def get_graph_number_cycles(cg:LigandNetwork, higher_bound:int=4)->int:
    graph = nx.DiGraph(cg.graph).to_undirected()
    raw_cycles = [c for c in nx.simple_cycles(graph, length_bound=higher_bound)]
    return len(raw_cycles)


def get_node_connectivities(cg: LigandNetwork)->dict[int, int]:
    return {n:sum([n in e for e in cg.edges]) for n in cg.nodes}


def get_norm_node_connectivities(cg: LigandNetwork)->dict[int, float]:
    n_edges = cg.graph.number_of_edges()
    return {n:sum([n in e for e in cg.edges])/n_edges for n in cg.nodes}


def get_node_scores(cg: LigandNetwork)->dict[int, float]:
    n_edges = cg.graph.number_of_edges()
    return {n:sum([e.annotations["score"] for e in cg.edges if(n in e)])/n_edges for n
            in cg.nodes}

def get_node_number_cycles(cg:LigandNetwork, higher_bound:int=4)->dict[int,
int]:
    graph = nx.DiGraph(cg.graph).to_undirected()
    raw_cycles = [n for c in nx.simple_cycles(graph, length_bound=higher_bound) for n in c]
    uni, cou = np.unique(raw_cycles, return_counts=True)
    return dict(zip(uni, cou))


def get_edge_failure_robustness(network, failure_rate=0.05, nrepeats=100):

    edges = list(network.edges)
    npics = int(np.round(len(network.edges) * failure_rate)) if (
            np.round(len(network.edges) * failure_rate) > 1) else 1

    connected=[]
    for _ in range(nrepeats):
        nn = network
        for _ in range(npics):
            i = random.randrange(0, len(edges), 1)
            nn = tools.delete_transformation(nn, edges[i])
        connected.append(get_is_connected(nn))

    r = np.mean(list(map(float, connected)))
    return r

