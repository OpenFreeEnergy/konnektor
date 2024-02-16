import numpy as np
from gufe import LigandNetwork

def get_node_connectivities(cg: LigandNetwork)->dict[int, int]:
    return {n:sum([n in e for e in cg.edges]) for n in cg.nodes}

def get_norm_node_connectivities(cg: LigandNetwork)->dict[int, float]:
    n_edges = cg.graph.number_of_edges()
    return {n:sum([n in e for e in cg.edges])/n_edges for n in cg.nodes}
