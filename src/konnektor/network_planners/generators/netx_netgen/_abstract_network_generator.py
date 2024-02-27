import abc
import networkx as nx

class _AbstractNetworkGenerator(abc.ABC):

    def __int__(self):
        pass

    def __call__(self, *args, **kwargs):
        return self.generate_network(*args, **kwargs)

    @abc.abstractmethod
    def generate_network(self, edges, weights)->nx.Graph:
        pass