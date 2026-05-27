# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import itertools
import warnings
from collections import Counter
from collections.abc import Iterable

from gufe import AtomMapper, Component, LigandNetwork

from .._map_scoring import _score_mappings
from ._abstract_network_generator import NetworkGenerator


class ExplicitNetworkGenerator(NetworkGenerator):
    def __init__(
        self,
        mappers: AtomMapper | list[AtomMapper],
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
        edges: Iterable[tuple[Component, Component]],
        nodes: Iterable[Component] | None = None,
    ) -> LigandNetwork:
        """
        Create a network with explicitly-defined edges and nodes.
        The network can be defined by specifying only edges, in which case the nodes are implicitly added.

        Parameters
        ----------
        edges: Iterable[Tuple[Component, Component]]
            Planned edges that will be connected with mappings and scores.
            Each Tuple represents one edge.

        nodes: Iterable[Component] | None
            A list of nodes to be included in the network.
            Optional, since the network can be defined by specifying only edges.
            This is useful for adding isolated (unconnected) nodes.

        Returns
        -------
        LigandNetwork


        Warns
        -----
        Warning
            Raises a warning if the network is not connected as a single network.
        """

        mappings = _score_mappings(
            possible_edges=edges,
            scorer=self.scorer,
            mappers=self.mappers,
            n_processes=self.n_processes,
            show_progress=self.progress,
        )

        network = LigandNetwork(edges=mappings, nodes=nodes)
        if not network.is_connected():
            warnings.warn("Generated network is not connected as a single network.")

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
        components : list[Component]
            ``Component``/s to place into the network.

        indices : list[tuple[int, int]]
            Edges to form between the ``Components``, represented as tuples of indices of the list of ``Component``/s.
            e.g. `[(3, 4), ...]` will create an edge between the 3rd and 4th molecules
            (remember that Python uses 0-based indexing)

        Returns
        -------
        LigandNetwork

        Raises
        ------
        IndexError
            Throws an error if the ``indices`` specified are not present in ``components``.j
        """
        edges = []

        for i, j in indices:
            try:
                edges.append((components[i], components[j]))
            except IndexError:
                raise IndexError(
                    f"Invalid ligand index. Requested ({i}, {j}) for iterable of length {len(components)}. Please choose values in range 0-{len(components) - 1}."
                )

        return self.generate_ligand_network(edges=edges, nodes=components)

    def generate_network_from_names(
        self,
        components: list[Component],
        names: list[tuple[str, str]],
    ) -> LigandNetwork:
        """
        Generate a :class:`.LigandNetwork` by specifying edges as tuples of names.

        Parameters
        ----------
        components : list[Component]
          ``Component``/s to place into the network.
        mapper: AtomMapper
          the atom mapper to use to construct edges
        names : list of tuples of names
          the edges to form where the values refer to names of the small molecules,
          eg ``[('benzene', 'toluene'), ...]`` will create an edge between the
          molecule with names 'benzene' and 'toluene'

        Returns
        -------
        LigandNetwork

        Raises
        ------
        KeyError
          If a name in ``names`` is not present in ``components``.
        ValueError
          If multiple molecules have the same name (molecule names must be unique)
        """
        nm2comp = {c.name: c for c in components}

        if len(nm2comp) < len(components):
            dupes = Counter(c.name for c in components)
            dupe_names = [k for k, v in dupes.items() if v > 1]
            raise ValueError(f"Duplicate names: {dupe_names}")

        edges = []
        for nameA, nameB in names:
            try:
                edges.append((nm2comp[nameA], nm2comp[nameB]))
            except KeyError:
                badnames = [nm for nm in itertools.chain.from_iterable(names) if nm not in nm2comp]
                available = [ligand.name for ligand in components]
                raise KeyError(f"Invalid name(s) requested {badnames}.  Available: {available}")

        return self.generate_ligand_network(edges=edges, nodes=components)
