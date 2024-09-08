# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import numpy as np
import networkx as nx
from matplotlib import pyplot as plt

from gufe import LigandNetwork
from . import color_gradient, OFE_COLORS


def get_node_connectivities(cg: nx.Graph) -> list[int]:
    """The connectivity of each node"""
    return [sum([n in e for e in cg.edges]) for n in cg.nodes]


def draw_ligand_network(
    network: LigandNetwork,
    title: str = "",
    ax: plt.Axes = None,
    node_size: int = 2050,
    edge_width: int = 3,
    fontsize: int = 18,
) -> plt.Figure:
    """visualize a LigandNetwork as matplotlib plot, indicating the graph topology.

    Parameters
    ----------
    network: LigandNetwork
        The network to be visualized
    title: str, optional
        plot title
    ax: plt.Axes, optional
        don't build new figure, if axes given, but add img to axes
    node_size: int, optional
        size of the node visualization. (default 2050)
    edge_width: int, optional
        widht of drawn edges. (default 3)
    fontsize: int, optional
        fontsize of labels. (default 18)

    Returns
    -------
    plt.Figure
        return the matplotlib figure
    """
    ligands = list(network.nodes)
    edge_map = {(m.componentA.name, m.componentB.name): m for m in network.edges}
    edges = list(sorted(edge_map.keys()))
    weights = [edge_map[k].annotations["score"] for k in edges]

    # g = network.graph
    g = nx.Graph()
    for n in ligands:
        g.add_node(n.name)
    g.add_weighted_edges_from(
        ebunch_to_add=[(e[0], e[1], w) for e, w in zip(edges, weights)]
    )

    pos = nx.spring_layout(g, weight=1)

    if ax is None:
        fig, ax = plt.subplots(figsize=[16, 9])
    else:
        fig = None

    connectivities = np.array(get_node_connectivities(network))
    mixins = (
        np.clip(
            connectivities / (sum(connectivities) / len(connectivities)),
            a_min=0,
            a_max=2,
        )
        / 2
    )

    cs = list(map(lambda x: color_gradient(mix=x, hex=False), mixins))

    nx.draw_networkx(
        g,
        pos=pos,
        with_labels=True,
        ax=ax,
        node_size=node_size,
        width=edge_width,
        node_color=cs,
        edge_color=OFE_COLORS[3],
        font_color=[1, 1, 1],
    )
    ax.set_title(title, fontsize=fontsize)  # +" #edges "+str(len(g.edges))

    return fig
