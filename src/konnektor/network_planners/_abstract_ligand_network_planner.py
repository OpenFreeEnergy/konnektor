import abc
import logging
import itertools
from typing import Iterable, Tuple

from gufe import SmallMoleculeComponent
from konnektor.utils  import LigandNetwork

log = logging.getLogger(__name__)
#log.setLevel(logging.WARNING)

class easyLigandNetworkPlanner(abc.ABC):
    def __init__(self, mapper, scorer, network_generator):
        self.mapper = mapper
        self.scorer =  scorer
        self.network_generator = network_generator

    def __call__(self, *args, **kwargs):
        return self.generate_ligand_network(*args, **kwargs)

    def _input_generate_all_possible_mappings(self, ligands)->Tuple:
        #Todo: this should actually become the function of max network.

        # First create a network with all the proposed mappings (scored)
        mapping_generator = itertools.chain.from_iterable(
            self.mapper.suggest_mappings(molA, molB)
            for molA, molB in itertools.combinations(ligands, 2))
        mappings = [mapping.with_annotations({'score': self.scorer(mapping), 'type': "core"})
                    for mapping in mapping_generator]

        #mapping_scores = [mapping.annotations['score'] for mapping in mappings]
        return ligands, mappings,


    def generate_ligand_network(self, ligands)->LigandNetwork:
        """Plan a Network which connects all ligands with minimal cost

        Parameters
        ----------
        ligands : Iterable[SmallMoleculeComponent]
        the ligands to include in the Network
        """
        nodes, edges, weights = self._input_generate_all_possible_mappings(ligands=ligands)
        selected_edges = self.network_generator(nodes=nodes, edges=edges, weights=weights)

        return LigandNetwork(edges=selected_edges, nodes=ligands)

