import networkx as nx
from matplotlib import pyplot as plt

ofe_colors =  [(49 / 256, 57 / 256, 77 / 256, 1),  # Badass Blue
                  (184 / 256, 87 / 256, 65 / 256, 1),  # Feeling spicy
                  (0, 147 / 256, 132 / 256, 1),  # Feeling sick
                  (217 / 256, 196 / 256, 177 / 256, 1),  # Beastlygrey
                  (102 / 256, 102 / 256, 102 / 256, 1),  # Sandy Sergio
                  (0 / 256, 47 / 256, 74 / 256, 1), ]  # otherBlue]


def draw_ligand_network(network, title="", ax=None, node_size=2050, edge_width=3):
    ligands = list(network.nodes)
    edge_map = {(m.componentA.name, m.componentB.name): m for m in network.edges}
    edges = list(sorted(edge_map.keys()))
    weights = [edge_map[k].annotations['score'] for k in edges]

    g = nx.Graph()
    [g.add_node(n.name) for n in ligands]
    g.add_weighted_edges_from(ebunch_to_add=[(e[0], e[1], w) for e,w in zip(edges,weights)])


    if(ax is None):
        fig = plt.figure(figsize=[16,9])
        ax = fig.gca()
    else:
        fig=None


    nx.draw_networkx(g, with_labels=True, ax=ax, node_size=node_size, width=edge_width,
                     node_color=ofe_colors[-1], edge_color=ofe_colors[3], font_color=[1,1,1])
    ax.set_title(title+" Network #edges "+str(len(g.edges)))

    return fig