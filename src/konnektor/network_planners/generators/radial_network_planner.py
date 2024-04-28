from typing import Iterable

from gufe import Component, LigandNetwork, AtomMapper, AtomMappingScorer
from konnektor.network_planners.generators.netx_netgen import RadialNetworkGenerator

from ._abstract_ligand_network_generator import LigandNetworkGenerator
from .maximal_network_planner import MaximalNetworkGenerator


class StarLigandNetworkGenerator(LigandNetworkGenerator):

    def __init__(self, mapper: AtomMapper, scorer: AtomMappingScorer,
                 nprocesses: int = 1, _initial_edge_lister: LigandNetworkGenerator = None):
        """
        The Star Ligand Network Planner or Radial Ligand Network Planner, set's one ligand into the center of a graph and connects all other ligands to it.

        Parameters
        ----------
        mapper : AtomMapper
            the atom mapper is required, to define the connection between two ligands.
        scorer : AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        nprocesses: int, optional
            number of processes that can be used for the network generation. (default: 1)
        _initial_edge_lister: LigandNetworkPlanner, optional
            this LigandNetworkPlanner is used to give the initial set of edges. For standard usage, the Maximal NetworPlanner is used.
            However in large scale approaches, it might be interesting to use the heuristicMaximalNetworkPlanner.. (default: MaximalNetworkPlanner)
        """
        if _initial_edge_lister is None:
            _initial_edge_lister = MaximalNetworkGenerator(mapper=mapper, scorer=scorer, nprocesses=nprocesses)

        super().__init__(mapper=mapper, scorer=scorer,
                         network_generator=RadialNetworkGenerator(),
                         nprocesses=nprocesses,
                         _initial_edge_lister=_initial_edge_lister)

    def generate_ligand_network(self, components: Iterable[Component],
                                central_component: Component = None) -> LigandNetwork:
        """
        generate a star map network for the given compounds. if a central component is defined,
        the planning stage is shortcutted to only connect the ligands to the central component.

        Parameters
        ----------
        components: Iterable[Component]
            the components to be used for the LigandNetwork

        central_component: Component, optional
            the central component can be given, in order to shortcut the calculations and enforce the central ligand.

        Returns
        -------
        LigandNetwork
            a star like network.
        """
        components = list(components)

        # If central Ligand not defined, get it from full graph.
        if (central_component is None):
            # Full Graph Construction
            initial_network = self._initial_edge_lister.generate_ligand_network(
                components=components)
            mappings = initial_network.edges

            # Translate Mappings to graphable:
            edge_map = {(components.index(m.componentA), components.index(m.componentB)): m for m in mappings}
            edges = list(sorted(edge_map.keys()))
            weights = [edge_map[k].annotations['score'] for k in edges]

            rg = self.network_generator.generate_network(edges=edges, weights=weights)
            selected_mappings = [edge_map[k] for k in rg.edges]

        else:  # Given central ligands: less effort. - Trivial Case
            if self.scorer is None:
                scorer = lambda x: -1
            else:
                scorer = self.scorer

            mapping_generators = [self.mapper.suggest_mappings(
                central_component, molA) for molA in components]
            selected_mappings = [mapping.with_annotations({'score': scorer(mapping)})
                                 for mapping_generator in mapping_generators for
                                 mapping in mapping_generator]

        return LigandNetwork(edges=selected_mappings, nodes=components)


RadialLigandNetworkPlanner = StarLigandNetworkGenerator
