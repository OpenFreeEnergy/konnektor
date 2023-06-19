import abc
from openfe.setup import LigandNetwork
import network_planner
import itertools
from typing import Iterable, List, Tuple

from gufe import SmallMoleculeComponent

class easyLigandNetworkPlanner(abc.ABC):
    def __init__(self, mapper, scorer, network_planner):
        self.mapper = mapper
        self.scorer =  scorer
        self.network_planner = network_planner

    @abc.abstractmethod
    def _input_generate_all_possible_mappings(self, ligands)->Tuple:
        pass

    def generate_ligand_network(self, ligands)->LigandNetwork:
        """Plan a Network which connects all ligands with minimal cost

        Parameters
        ----------
        ligands : Iterable[SmallMoleculeComponent]
        the ligands to include in the Network
        """
        nodes, edges, weights = self._input_generate_all_possible_mappings(ligands=ligands)
        selected_edges = self.network_planner(nodes=nodes, edges=edges, weights=weights)

        return LigandNetwork(edges=selected_edges, nodes=ligands)

class MinimalSpanningTreeLigandNetworkPlanner(easyLigandNetworkPlanner):

    def __init__(self, mapper, scorer):
        """Plan a Network which connects all ligands with minimal cost

        Parameters
        ----------
        mappers : Iterable[AtomMapper]
        the AtomMappers to use to propose mappings.  At least 1 required,
        but many can be given, in which case all will be tried to find the
        lowest score edges
        scorer : Scoring function
        any callable which takes a AtomMapping and returns a float
        """
        super.__init__(mapper=mapper, scorer=scorer,
                       network_planner=network_planner.mst_network_planner())

    def _input_generate_all_possible_mappings(self, ligands)->Tuple:
        # First create a network with all the proposed mappings (scored)
        mapping_generator = itertools.chain.from_iterable(
            self.mapper.suggest_mappings(molA, molB)
            for molA, molB in itertools.combinations(ligands, 2)
        )
        mappings = [mapping.with_annotations({'score': self.scorer(mapping)})
                    for mapping in mapping_generator]

        mapping_scores = [mapping.score for mapping in mappings]
        return ligands, mappings,

    def generate_ligand_network(self, ligands) ->LigandNetwork:

        mappings = self._input_generate_all_possible_mappings(ligands=ligands)
        index_nodes = {ligands.index(l): l for l in ligands}
        edge_map = {(ligands.index(m.componentA), ligands.index(m.componentA)): m for m in mappings}
        edges = list(sorted(edge_map.keys()))
        weights = [edges[k].score for k in sorted(edge_map)]

        mg = self.network_planner.generate_network(edges, weights)

        if(len(mg.nodes) < len(ligands)):
            missing_nodes = [l for l in ligands if(index_nodes[l] in mg.nodes)]
            raise RuntimeError("Unable to create edges to some nodes: "
                               + str(list(missing_nodes)))

        selected_mappings = [edge_map[k] for k in mg.edges]

        return LigandNetwork(edges=selected_mappings, nodes=ligands)


class CyclicLigandNetworkPlanner(easyLigandNetworkPlanner):

    def __init__(self, mapper, scorer, node_present_in_cycles=2, cycle_sizes:Union[int, List[int]]=3 ):
        super.__init__(mapper=mapper, scorer=scorer, network_planner=network_planner.cyclic_network_planner(node_cycle_connectivity=node_present_in_cycles, sub_cycle_size_range=cycle_sizes))

    def __call__(self, *args, **kwargs):
        self.generate_ligand_network(*args, **kwargs)

    def generate_ligand_network(self, ligands: Iterable[SmallMoleculeComponent]) ->LigandNetwork:
            # Build Full Graph

            ligands = list(ligands)

            mapping_generator = itertools.chain.from_iterable(
                self.mapper.suggest_mappings(molA, molB)
                for molA, molB in itertools.combinations(ligands, 2))
            mappings = [mapping.with_annotations({'score': self.scorer(mapping)})
                        for mapping in mapping_generator]

            # get cyclic_graph
            edge_map = {(ligands.index(m.componentA), ligands.index(m.componentA)): m for m in mappings}
            edges = list(sorted(edge_map.keys()))
            weights = [edges[k].score for k in sorted(edge_map)]

            cg = self.network_planner.generate_network(edges=edges, weights=weights)

            selected_mappings = [edge_map[k] for k in cg.edges]
            return LigandNetwork(edges=selected_mappings, nodes=ligands)


class RadialLigandNetworkPlanner(easyLigandNetworkPlanner):

    def __init__(self, mapper, scorer):
        super.__init__(mapper=mapper, scorer=scorer, network_planner=network_planner.radial_network_planner())

    def __call__(self, *args, **kwargs):
        self.generate_ligand_network(*args, **kwargs)

    def generate_ligand_network(self, ligands: Iterable[SmallMoleculeComponent], central_ligand=None) ->LigandNetwork:
            # Build Full Graph
            ligands = list(ligands)

            if(central_ligand is None): #Full Graph
                mapping_generator = itertools.chain.from_iterable(
                    self.mapper.suggest_mappings(molA, molB)
                    for molA, molB in itertools.combinations(ligands, 2))
                mappings = [mapping.with_annotations({'score': self.scorer(mapping)})
                            for mapping in mapping_generator]

                #Translate Mappings to graphable:
                edge_map = {(ligands.index(m.componentA), ligands.index(m.componentA)): m for m in mappings}
                edges = list(sorted(edge_map.keys()))
                weights = [edges[k].score for k in sorted(edge_map)]

                rg = self.network_planner.generate_network(edges=edges, weights=weights)

                selected_mappings = [edge_map[k] for k in rg.edges]

            else:   #Given central ligands: less effort. - Trivial Case
                mapping_generator = [self.mapper.suggest_mappings(central_ligand, molA) for molA in ligands]
                selected_mappings = [mapping.with_annotations({'score': self.scorer(mapping)})
                            for mapping in mapping_generator]

            return LigandNetwork(edges=selected_mappings, nodes=ligands)


'''
class ligandNetworkPlanner():
    def __int__(self, mappers=[], scorers=[], network_planner=None):
        self.mappers = mappers
        self.mapping_scorer = scorers
        self.network_planner = network_planner

    def __call__(self, *args, **kwargs):
        self.generate_ligand_network()

    def generate_ligand_network(self, ligands:Iterable):

        mappings = [[[mapper(molA, molB) for molB in ligands[i:]] for i, molA in enumerate(ligands)] for mapper in self.mappers]
        scores = [[[(mapping_scorer(mapping) for mapping_scorer in self.mapping_scorer) for mapping in molA_mappings] for molA_mappings in mapper_mappings] for mapper_mappings in mappings]

        #flatten_stuff
        mappings_scores =[]
        for mapping_mappper, mapping_scores in zip(mappings, scores):
            for mol_mappings, mol_score in zip(mapping_mappper, mapping_scores):
                mappings_scores.extend(list(zimp(mol_mappings, mol_score)))

        self.network_planner()

'''