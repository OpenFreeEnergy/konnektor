# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import networkx as nx
import numpy as np
from urllib import parse

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from ipywidgets import interact, fixed
import ipycytoscape

import gufe
from gufe import AtomMapping, LigandNetwork
from gufe.visualization.mapping_visualization import draw_mapping

from . import color_gradient


def get_node_connectivities(cg: LigandNetwork) -> list[int]:
    """The connectivity of each node"""
    return [sum([n in e for e in cg.edges]) for n in cg.nodes]


# some code borrowed from pen:
# https://iwatobipen.wordpress.com/2020/03/30/draw-scaffold-tree-as-network-with-molecular-image-rdkit-cytoscape/
def mol2svg(mol: Chem.Mol) -> str:
    try:
        Chem.rdmolops.Kekulize(mol)
    except:
        pass
    drawer = rdMolDraw2D.MolDraw2DSVG(350, 300)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer,
                                       mol)  # , legend=mol.GetProp("_Name"))
    drawer.SetColour(
        (184 / 256, 87 / 256, 65 / 256))  # Transparent white background
    drawer.FinishDrawing()

    svg = drawer.GetDrawingText()
    impath = 'data:image/svg+xml;charset=utf-8,' + parse.quote(svg, safe="")

    return impath


def map2svg(mapping: AtomMapping) -> str:
    grid_x, grid_y = {1: (1, 1), 2: (2, 1), }[2]
    d2d = rdMolDraw2D.MolDraw2DSVG(
        grid_x * 300, grid_y * 300, 300, 300)
    d2d.SetColour((1, 1, 1))

    svg = draw_mapping(mol1=mapping.componentA.to_rdkit(),
                       mol2=mapping.componentB.to_rdkit(),
                       mol1_to_mol2=mapping.componentA_to_componentB,
                       d2d=d2d)
    impath = 'data:image/svg+xml;charset=utf-8,' + parse.quote(svg, safe="")

    return impath


def _build_cytoscape(network: gufe.LigandNetwork, layout: str = "concentric",
                     show_molecules: bool = True,
                     show_mappings: bool = False) -> ipycytoscape.CytoscapeWidget:
    ligands = list(network.nodes)
    edge_map = {(m.componentA.name, m.componentB.name): m for m in
                network.edges}
    edges = list(sorted(edge_map.keys()))
    weights = [edge_map[k].annotations['score'] for k in edges]

    connectivities = np.array(get_node_connectivities(network))
    if not len(connectivities):
        mixins = np.array([0])
        cs = list(map(lambda x: color_gradient(mix=x), mixins))
    else:
        mixins = np.clip(
            connectivities / (sum(connectivities) / len(connectivities)),
            a_min=0, a_max=2) / 2

        cs = list(map(lambda x: color_gradient(mix=x), mixins))

    g = nx.Graph()

    g.add_nodes_from(
        (n.name,
         {'name': n.name, 'classes': 'ligand',
          'img': mol2svg(n.to_rdkit()), 'col': c},
         )
        for n, c in zip(ligands, cs)
    )

    if show_mappings:
        g.add_nodes_from(
            (f'{e[0]}-{e[1]}',
             {'name': f'{e[0]}-{e[1]}', 'classes': 'mapping',
              'lab': f'{e[0]} - {e[1]}\nscore: {w:2.2F}',
              'weight': f'{w:2.2F}', 'img': map2svg(edge_map[e])},
             )
            for e, w in zip(edges, weights)
        )
    else:
        g.add_nodes_from(
            (f'{e[0]}-{e[1]}',
             {'name': f'{e[0]}-{e[1]}', 'classes': 'mapping',
              'lab': f'{e[0]} - {e[1]}\nscore: {w:2.2F}',
              'weight': f'{w:2.2F}',},
             )
            for e, w in zip(edges, weights)
        )

    for e, w in zip(edges, weights):
        g.add_edge(e[0],f"{e[0]}-{e[1]}")
        g.add_edge(f"{e[0]}-{e[1]}", e[1])

    g = g.to_undirected()

    # Styling
    node_style_ligand = {
        'selector': 'node.ligand',
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
    if show_molecules:
        node_style_ligand["css"]["background-image"] = "data(img)"
        node_style_ligand["css"]["text-valign"] = "top"
        node_style_ligand["css"]["width"] = 350
        node_style_ligand["css"]["height"] = 300

    node_style_mapping = {
        'selector': 'node.mapping',
        'css': {
            'content': 'data(lab)',
            'text-valign': 'center',
            'text-outline-width': 1,
            'text-outline-color': '#31394d',
            "text-wrap": "wrap",
            "font-size": 20,
            'color': '#31394d',
            'shape': 'roundrectangle',
            'width': 150,
            'height': 100,
            'nested': True,
            'background-color': '#d9c4b1',
            'background-fit': 'contain',
            'border-width': 5,
            'border-color': '#d9c4b1',
        }
    }
    if show_molecules and show_mappings:
        node_style_mapping["css"]["background-image"] = "data(img)"
        node_style_mapping["css"]["background-color"] = "white"
        node_style_mapping["css"]["text-valign"] = "top"
        node_style_mapping["css"]["width"] = 300
        node_style_mapping["css"]["height"] = 150

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
            'tooltip': 'data(weight)',
            'tooltip-background-image': 'data(edge_img)',
        }
    }

    undirected = ipycytoscape.CytoscapeWidget(g, )
    undirected.wheel_sensitivity = 0.2
    undirected.set_tooltip_source("name")
    undirected.set_style([node_style_ligand, node_style_mapping, edge_style])
    undirected.set_layout(name=layout, nodeSpacing=50, edgeLengthVal=50)
    return undirected


def draw_network_widget(network: gufe.LigandNetwork, layout: str = "cose",
                        show_molecules: bool = True,
                        show_mappings: bool = False) -> ipycytoscape.CytoscapeWidget:
    """For use in a jupyter noterbook, visualise a LigandNetwork

    Parameters
    ----------
    network: gufe.LigandNetwork
      the network to visualise
    layout : str, optional
      how to initially layout the nodes, can be one of X/Y/Z
      defaults to 'cose'
    show_molecule: bool, optional
      if to show molecule images on the representation, default True
    show_mappings: bool, optional
      if to show mapping images on the representation, default False
    """

    @interact(network=fixed(network), layout=['dagre', 'cola', 'breadthfirst',
                                              'concentric', 'cose'])
    def interactive_widget(network=network, layout=layout,
                           show_molecules=show_molecules,
                           show_mappings=show_mappings):
        v = _build_cytoscape(network=network, layout=layout,
                             show_molecules=show_molecules,
                             show_mappings=show_mappings)
        return v

    return interactive_widget(network=network, layout=layout,
                              show_molecules=show_molecules)
