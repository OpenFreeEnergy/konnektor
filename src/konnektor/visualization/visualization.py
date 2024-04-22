import numpy as np
import networkx as nx
from matplotlib import pyplot as plt

ofe_colors =  [(49 / 256, 57 / 256, 77 / 256),  # Badass Blue
              (184 / 256, 87 / 256, 65 / 256),  # Feeling spicy
              (0, 147 / 256, 132 / 256),  # Feeling sick
              (217 / 256, 196 / 256, 177 / 256),  # Beastlygrey
              (217 / 256, 196 / 256, 177 / 256),  # Sandy Sergio
               (238/256, 192/256, 68/256), #GOld
              (0 / 256, 47 / 256, 74 / 256), ]  # otherBlue]

def color_gradient(c1=ofe_colors[1], c2=ofe_colors[2], c3=ofe_colors[1], mix=0):
    c1=np.array(c1)
    c2=np.array(c2)
    c3=np.array(c3)
    mix = np.array(mix, ndmin=1)

    if(mix > 0.5):
        m = mix-0.5
        c = (0.5-m)*c2 + m*c3
    else:
        m = mix
        c = (0.5-m)*c1 + m*c2
    return c

def get_node_connectivities(cg):
    return [sum([n in e for e in cg.edges]) for n in cg.nodes]



def draw_ligand_network(network, title="", ax=None, node_size=2050, edge_width=3, fontsize=18):
    ligands = list(network.nodes)
    edge_map = {(m.componentA.name, m.componentB.name): m for m in network.edges}
    edges = list(sorted(edge_map.keys()))
    weights = [edge_map[k].annotations['score'] for k in edges]

    #g = network.graph
    g = nx.Graph()
    [g.add_node(n.name) for n in ligands]
    g.add_weighted_edges_from(ebunch_to_add=[(e[0], e[1], w) for e,w in zip(
     edges,weights)])

    if(ax is None):
        fig, ax = plt.subplots(figsize=[16,9])
    else:
        fig=None

    connectivities = np.array(get_node_connectivities(network))
    mixins = np.clip(connectivities / (sum(connectivities)/len(connectivities)), a_min=0, a_max=2)/2

    cs = list(map(lambda x: color_gradient(mix=x), mixins))

    nx.draw_networkx(g, with_labels=True, ax=ax, node_size=node_size, width=edge_width,
                     node_color=cs, edge_color=ofe_colors[3], font_color=[1,1,1])
    ax.set_title(title, fontsize=fontsize) #+" #edges "+str(len(g.edges))

    return fig
