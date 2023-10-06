import abc

from gufe import network

class Network():
    pass

class _AbstractNetworkGenerator(abc.ABC):

    def __int__(self):
        pass

    def __call__(self, *args, **kwargs):
        return self.generate_network(*args, **kwargs)

    @abc.abstractmethod
    def generate_network(self, edges, weights)->Network:
        pass