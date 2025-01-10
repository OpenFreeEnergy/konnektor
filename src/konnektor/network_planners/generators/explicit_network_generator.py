# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import warnings
import itertools
from collections import Counter
from typing import Iterable, Tuple, Union

from gufe import Component, LigandNetwork, AtomMapper

from ._abstract_network_generator import NetworkGenerator
from ._parallel_mapping_pattern import _parallel_map_scoring


class ExplicitNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: Union[AtomMapper, list[AtomMapper]],
        scorer,
        n_processes: int = 1,
        progress: bool = False,
    ):
        """

        Parameters
        ----------
        mapper: AtomMapper
            Defines the connection between two ligands.
        scorer: AtomMappingScorer
            scoring function evaluating an atom mapping, and giving a score between [0,1].
        n_processes: int
            number of processes to use to build the ligand network
        progress: bool, optional
            if true a progress bar will be displayed. (default: False)
        """
        super().__init__(
            mappers=mappers,
            scorer=scorer,
            n_processes=n_processes,
            progress=progress,
            network_generator=None,
        )

    def generate_ligand_network(
        self,
        edges: Iterable[Tuple[Component, Component]],
    ) -> LigandNetwork:
        """
        Create a network with pre-defined edges.

        This can be used as initial_edge_lister

        Parameters
        ----------
        edges: Iterable[Tuple[Component,Component]]
            planned edges, that will be connected with mappings and scores.
            Each Tuple in this case represent one edge.

        Returns
        -------
        LigandNetwork
            the provided network.
        """
        nodes = list(set([n for e in edges for n in e]))

        mappings = _parallel_map_scoring(
            possible_edges=edges,
            scorer=self.scorer,
            mappers=self.mappers,
            n_processes=self.n_processes,
            show_progress=self.progress,
        )

        network = LigandNetwork(edges=mappings, nodes=nodes)
        if not network.is_connected():
            warnings.warn("Generated network is not fully connected")

        return network

    def generate_network_from_indices(
        self,
        components: list[Component],
        indices: list[tuple[int, int]],
    ) -> LigandNetwork:
        """
        Generate a :class:`.LigandNetwork` by specifying edges as tuples of indices.

        Parameters
        ----------
        components : list of Components
          the small molecules to place into the network
        mapper: AtomMapper
          the atom mapper to use to construct edges
        indices : list of tuples of indices
          the edges to form where the values refer to names of the small molecules,
          eg `[(3, 4), ...]` will create an edge between the 3rd and 4th molecules
          remembering that Python uses 0-based indexing

        Returns
        -------
        LigandNetwork

        Raises
        ------
        IndexError
          if an invalid ligand index is requested
        """
        edges = []
        for i, j in indices:
            try:
                edges.append((components[i], components[j]))
            except IndexError:
                raise IndexError(
                    f"Invalid ligand id, requested {i} {j} "
                    f"with {len(components)} available"
                )

        return self.generate_ligand_network(edges=edges)

    def generate_network_from_names(
        self,
        components: list[Component],
        names: list[tuple[str, str]],
    ) -> LigandNetwork:
        """
        Generate a :class:`.LigandNetwork` by specifying edges as tuples of names.

        Parameters
        ----------
        components : list of Components
          the small molecules to place into the network
        mapper: AtomMapper
          the atom mapper to use to construct edges
        names : list of tuples of names
          the edges to form where the values refer to names of the small molecules,
          eg `[('benzene', 'toluene'), ...]` will create an edge between the
          molecule with names 'benzene' and 'toluene'

        Returns
        -------
        LigandNetwork

        Raises
        ------
        KeyError
          if an invalid name is requested
        ValueError
          if multiple molecules have the same name (this would otherwise be
          problematic)
        """
        nm2comp = {l.name: l for l in components}

        if len(nm2comp) < len(components):
            dupes = Counter((l.name for l in components))
            dupe_names = [k for k, v in dupes.items() if v > 1]
            raise ValueError(f"Duplicate names: {dupe_names}")

        edges = []
        for nameA, nameB in names:
            try:
                edges.append((nm2comp[nameA], nm2comp[nameB]))
            except KeyError:
                badnames = [
                    nm
                    for nm in itertools.chain.from_iterable(names)
                    if nm not in nm2comp
                ]
                available = [ligand.name for ligand in components]
                raise KeyError(
                    f"Invalid name(s) requested {badnames}.  Available: {available}"
                )

        return self.generate_ligand_network(edges=edges)
