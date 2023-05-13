import abc

from gufe import network

class Network():
    pass

class _AbstractNetworkPlanner(abc.ABC):

    def __int__(self):
        pass

    def __call__(self, *args, **kwargs):
        self.generate_network(*args, **kwargs)

    @abc.abstractmethod
    def generate_network(self)->Network:
        pass