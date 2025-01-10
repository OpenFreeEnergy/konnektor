# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor

import abc

import networkx as nx


class _AbstractNetworkAlgorithm(abc.ABC):
    def __call__(self, *args, **kwargs):
        return self.generate_network(*args, **kwargs)

    @abc.abstractmethod
    def generate_network(self, edges: list[tuple[int, int]], weights: list[float]) -> nx.Graph:
        raise NotImplementedError()


class _AbstractNetworkConcatenator(abc.ABC):
    def __call__(self, *args, **kwargs):
        return self.concatenate_networks(*args, **kwargs)

    @abc.abstractmethod
    def concatenate_networks(
        self,
        nodesA: list[int],
        nodesB: list[int],
        edges: list[tuple[int, int]],
        weights: list[float],
    ) -> nx.Graph:
        raise NotImplementedError()
