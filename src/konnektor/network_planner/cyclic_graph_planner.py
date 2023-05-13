import itertools

from typing import Iterable, Callable, Union, Optional

from gufe import SmallMoleculeComponent, AtomMapper

from openfe.setup.ligand_network import LigandNetwork    # only temproary
from ._abstract_network_planner import _AbstractNetworkPlanner, Network


class cyclicGraph(_AbstractNetworkPlanner):

    def __init__(self, node_connectivity:int=2,
                 sub_cycle_size_range:Iterable[int]=[3]):
        self.node_connectivity=node_connectivity
        self.sub_cycle_size_range = sub_cycle_size_range


    def generate_network(self) ->Network:
