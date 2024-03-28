import networkx as nx
import numpy as np
from urllib import parse

from ipywidgets import interact, fixed
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import ipycytoscape

from gufe.visualization.mapping_visualization import draw_mapping

ofe_colors = [(49 / 256, 57 / 256, 77 / 256),  # Badass Blue
              (184 / 256, 87 / 256, 65 / 256),  # Feeling spicy
              (0, 147 / 256, 132 / 256),  # Feeling sick
              (217 / 256, 196 / 256, 177 / 256),  # Beastlygrey
              (217 / 256, 196 / 256, 177 / 256),  # Sandy Sergio
              (238 / 256, 192 / 256, 68 / 256),  # GOld
              (0 / 256, 47 / 256, 74 / 256), ]  # otherBlue]

rgb2hex = lambda r, g, b: '#%02x%02x%02x' % (int(r * 256), int(g * 256), int(b * 256))


def color_gradient(c1=ofe_colors[1], c2=ofe_colors[2], c3=ofe_colors[1], mix=0):
    c1 = np.array(c1)
    c2 = np.array(c2)
    c3 = np.array(c3)
    mix = np.array(mix, ndmin=1)

    if (mix > 0.5):
        m = mix - 0.5
        c = (0.5 - m) * c2 + m * c3
    else:
        m = mix
        c = (0.5 - m) * c1 + m * c2
    return rgb2hex(*c)


def get_node_connectivities(cg):
    return [sum([n in e for e in cg.edges]) for n in cg.nodes]


# some code borrowed from pen:
# https://iwatobipen.wordpress.com/2020/03/30/draw-scaffold-tree-as-network-with-molecular-image-rdkit-cytoscape/
def mol2svg(mol):
    try:
        Chem.rdmolops.Kekulize(mol)
    except:
        pass
    drawer = rdMolDraw2D.MolDraw2DSVG(350, 300)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)  # , legend=mol.GetProp("_Name"))
    drawer.SetColour((184 / 256, 87 / 256, 65 / 256))  # Transparent white background
    drawer.FinishDrawing()

    svg = drawer.GetDrawingText()
    impath = 'data:image/svg+xml;charset=utf-8,' + parse.quote(svg, safe="")

    return impath


def build_cytoscape(network, layout="concentric", show_mol_img=True):
    ligands = list(network.nodes)
    edge_map = {(m.componentA.name, m.componentB.name): m for m in network.edges}
    edges = list(sorted(edge_map.keys()))
    weights = [edge_map[k].annotations['score'] for k in edges]
    imgs = [edge_map[k] for k in edges]
    connectivities = np.array(get_node_connectivities(network))
    mixins = np.clip(connectivities / (sum(connectivities) / len(connectivities)), a_min=0, a_max=2) / 2

    cs = list(map(lambda x: color_gradient(mix=x), mixins))

    # build a graph
    g = nx.Graph()
    [g.add_node(n.name, name=n.name, img=mol2svg(n.to_rdkit()), col=c) for n, c in zip(ligands, cs)]
    g.add_weighted_edges_from(
        ebunch_to_add=[(e[0], e[1], "{:2.2F}".format(w),) for e, w, i in zip(edges, weights, imgs)])
    g = g.to_undirected()

    # Styling
    node_style = {
        'selector': 'node',
        'css': {
            'content': 'data(name)',
            'text-valign': 'center',
            'text-outline-width': 2,
            'text-outline-color': '#31394d',
            "font-size": 32,
            'color': '#d9c4b1',
            'shape': 'circle',
            'width': 350,
            'height': 300,
            'background-color': 'data(col)',
            'background-fit': 'contain',
            'border-width': 10,
            'border-color': 'data(col)',
        }
    }
    if show_mol_img:
        node_style["css"]["background-image"] = "data(img)"
        node_style["css"]["text-valign"] = "top"
        node_style["css"]["width"] = 350
        node_style["css"]["height"] = 300

    edge_style = {
        'selector': 'edge',
        'style': {
            'width': 4,
            'line-color': '#d9c4b1',
            'curve-style': 'bezier',
            'label': 'data(weight)',
            'color': '#31394d',
            "font-size": 26,
            'text-background-color': '#d9c4b1',  # Set the background color
            'text-background-opacity': 0.7,  # Adjust opacity (0 to 1)
        }
    }

    selected_style = {
        'selector': ':selected',
        'css': {
            'background-color': '#ede3da',
            'line-color': '#ede3da',
            'target-arrow-color': '#ede3da',
            'source-arrow-color': '#ede3da',
            'text-outline-color': '#ede3da'
        }
    }

    undirected = ipycytoscape.CytoscapeWidget()
    undirected.graph.add_graph_from_networkx(g)
    undirected.set_style([node_style, edge_style, selected_style])
    undirected.set_layout(name=layout, nodeSpacing=50, edgeLengthVal=50)
    return undirected


def draw_network_widget(network, layout="concentric", show_mol_img=True):
    @interact(network=fixed(network), layout=['dagre', 'breadthfirst', 'concentric', 'cose'],
              show_mol_img=[True, False])
    def interactive_widget(network=network, layout=layout, show_mol_img=show_mol_img):
        v = build_cytoscape(network=network, layout=layout, show_mol_img=show_mol_img)
        return v

    return interactive_widget(network=network, layout=layout, show_mol_img=show_mol_img)


