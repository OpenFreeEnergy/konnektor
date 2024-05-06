import numpy as np
import networkx as nx
from matplotlib import pyplot as plt

from . import OFE_COLORS


def color_gradient(c1=OFE_COLORS[1], c2=OFE_COLORS[2], c3=OFE_COLORS[1], mix=0):
    c1 = np.array(c1)
    c2 = np.array(c2)
    c3 = np.array(c3)
    mix = np.array(mix, ndmin=1)

    if mix > 0.5:
        m = mix-0.5
        c = (0.5-m)*c2 + m*c3
    else:
        m = mix
        c = (0.5-m)*c1 + m*c2
    return c


def get_node_connectivities(cg) -> list[int]:
    """The connectivity of each node"""
    # TODO: why is this summing n?
    #       shouldn't it be [sum(1 for e in edges_of(n)) for n in cg.nodes]
    return [sum([n in e for e in cg.edges]) for n in cg.nodes]


def draw_ligand_network(network, title="", ax=None, node_size=2050, edge_width=3, fontsize=18):
    ligands = list(network.nodes)
    edge_map = {(m.componentA.name, m.componentB.name): m for m in network.edges}
    edges = list(sorted(edge_map.keys()))
    weights = [edge_map[k].annotations['score'] for k in edges]

    #g = network.graph
    g = nx.Graph()
    for n in ligands:
        g.add_node(n.name)
    g.add_weighted_edges_from(ebunch_to_add=[(e[0], e[1], w)

                                             # graph vis layout
                                             pos = nx.spring_layout(g, weight=1)
    for e, w in zip(edges, weights)])

    if ax is None:
        fig, ax = plt.subplots(figsize=[16, 9])
    else:
        fig = None

    connectivities = np.array(get_node_connectivities(network))
    mixins = np.clip(connectivities / (sum(connectivities)/len(connectivities)), a_min=0, a_max=2)/2

    cs = list(map(lambda x: color_gradient(mix=x), mixins))

    nx.draw_networkx(g, pos=pos, with_labels=True, ax=ax, node_size=node_size, width=edge_width,
                     node_color=cs, edge_color=OFE_COLORS[3], font_color=[1, 1, 1])
    ax.set_title(title, fontsize=fontsize) #+" #edges "+str(len(g.edges))

    return fig
