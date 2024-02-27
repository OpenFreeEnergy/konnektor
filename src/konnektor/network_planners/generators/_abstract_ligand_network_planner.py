import abc
import logging
from typing import Iterable

from gufe import SmallMoleculeComponent, LigandNetwork

log = logging.getLogger(__name__)
#log.setLevel(logging.WARNING)


class LigandNetworkPlanner(abc.ABC):
    progress: bool = False
    nprocesses: int

    def __init__(self, mapper, scorer, network_generator, nprocesses:int=1,
                 _initial_edge_lister=None):
        self.mapper = mapper
        self.scorer =  scorer
        self.network_generator = network_generator
        self.nprocesses=nprocesses
        self._initial_edge_lister = _initial_edge_lister

    def __call__(self, *args, **kwargs)-> LigandNetwork:
        return self.generate_ligand_network(*args, **kwargs)

    def generate_ligand_network(self, ligands)->LigandNetwork:
        """Plan a Network which connects all ligands with minimal cost

        Parameters
        ----------
        ligands : Iterable[SmallMoleculeComponent]
        the ligands to include in the Network
        """
        initial_network = self._initial_edge_lister.generate_ligand_network(
            nodes=ligands)

        nodes = ligands
        edge_map = {(ligands.index(m.componentA), ligands.index(
            m.componentB)): m for m in initial_network.edges}
        edges = list(edge_map.keys())
        weights = [edge_map[k].annotations['score'] for k in edges]

        selected_edges = self.network_generator(nodes=nodes, edges=edges, weights=weights)

        return LigandNetwork(edges=selected_edges, nodes=ligands)

