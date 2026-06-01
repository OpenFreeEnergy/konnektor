# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/konnektor


import warnings

import networkx as nx
import numpy as np
from gufe import LigandNetwork, SmallMoleculeComponent

from .. import network_tools as tools


#  TODO: deprecate this?
def get_is_connected(ligand_network: LigandNetwork) -> bool:
    """
    Check whether all nodes in the LigandNetwork are connected to each other, ignoring edge direction.

    A False value indicates that either some nodes have no edges
    or that there are separate networks that do not link to each other.

    .. note ::
        This function is deprecated, please use the ``.is_connected`` method of the ``ligand_network`` instead.
        See `gufe.LigandNetwork <https://gufe.openfree.energy/en/latest/generated/gufe.ligandnetwork.html#gufe.ligandnetwork.LigandNetwork.is_connected>`_ docs for more details.

    Parameters
    ----------
    ligand_network: LigandNetwork

    Returns
    -------
    bool
        True if the LigandNetwork is connected, otherwise False.
    """
    msg = "`get_is_connected` is deprecated, please use the ligand_network's `.is_connected()` method instead"
    warnings.warn(msg, DeprecationWarning)

    return ligand_network.is_connected()


def get_network_score(ligand_network: LigandNetwork) -> float:
    """
    Calculate the network score based on summation of the edge scores.

    Parameters
    ----------
    ligand_network: LigandNetwork

    Returns
    -------
    float
        Summation of all edge scores in the LigandNetwork.
    """
    # TODO: "score" is not guaranteed to exist as an annotation,
    # so we should have better error handling here.
    score = sum([e.annotations["score"] for e in ligand_network.edges])
    return score


def get_network_cost(ligand_network: LigandNetwork) -> float:
    """
    Calculate the network cost based on summation of the edge costs,
    `cost = 1 - score`.

    Parameters
    ----------
    ligand_network: LigandNetwork

    Returns
    -------
    float
        Summation of all edge costs in the LigandNetwork.
    """
    score = sum([1 - float(e.annotations["score"]) for e in ligand_network.edges])
    return score


def get_network_efficiency(ligand_network: LigandNetwork) -> float:
    """
    Calculate the network efficiency, defined as the network score normalized by the number of edges in the network.

    Parameters
    ----------
    ligand_network: LigandNetwork

    Returns
    -------
    float
        Network efficiency score.
    """
    score = sum([e.annotations["score"] for e in ligand_network.edges]) / len(ligand_network.edges)
    return score


def get_number_of_network_cycles(ligand_network: LigandNetwork, higher_bound: int = 3) -> int:
    """Calculate the number of cycles present in the network for cycles up to size `higher_bound`.

    Parameters
    ----------
    ligand_network : LigandNetwork
    higher_bound : int, optional
        Largest number of nodes in cycle to count, by default 3.

    Returns
    -------
    int
        Number of cycles of size <= `higher_bound` present in the LigandNetwork.
    """
    graph = nx.DiGraph(ligand_network.graph).to_undirected()
    raw_cycles = [str(sorted(c)) for c in nx.simple_cycles(graph, length_bound=higher_bound)]
    return len(raw_cycles)


def get_component_connectivities(
    ligand_network: LigandNetwork,
    normalize: bool = False,
) -> dict[SmallMoleculeComponent, float | int]:
    """Calculate the connectivity for all nodes in the graph.

    Parameters
    ----------
    ligand_network : LigandNetwork
    normalize : bool, optional
        If True, normalize each node's connectivity by the total number of edges in the LigandNetwork, by default False.

    Returns
    -------
    dict[SmallMoleculeComponent, float | int]
        Dict of nodes and their corresponding connectivity (`int` if `normalize=False`, otherwise `float`).
    """
    if normalize:
        n_edges = ligand_network.graph.number_of_edges()
        return {
            n: sum([n in e for e in ligand_network.edges]) / n_edges for n in ligand_network.nodes
        }
    else:
        return {n: sum([n in e for e in ligand_network.edges]) for n in ligand_network.nodes}


def get_component_scores(
    ligand_network: LigandNetwork,
    normalize: bool = True,  # TODO: why is this True by default?
) -> dict[SmallMoleculeComponent, float]:
    """Calculate the score for all nodes in the graph, where a node's score is defined as
    the sum of the scores of all its edges.

    Parameters
    ----------
    ligand_network : LigandNetwork

    normalize : bool, optional
        If True, normalize each node's score by the total number of edges in the LigandNetwork, by default True.

    Returns
    -------
    dict[SmallMoleculeComponent, float]
        Dict of nodes (SmallMoleculeComponent) and their corresponding scores.
    """
    if normalize:
        n_edges = ligand_network.graph.number_of_edges()
        return {
            n: sum([e.annotations["score"] for e in ligand_network.edges if n in e]) / n_edges
            for n in ligand_network.nodes
        }
    else:
        return {
            n: sum([e.annotations["score"] for e in ligand_network.edges if n in e])
            for n in ligand_network.nodes
        }


def get_component_number_cycles(
    ligand_network: LigandNetwork,
    higher_bound: int = 4,
) -> dict[SmallMoleculeComponent, int | np.int64]:
    """Calculate the number of cycles each node is present in.

    Parameters
    ----------
    ligand_network : LigandNetwork
    higher_bound : int, optional
        Largest cycle (by number of nodes) to be considered a "cycle" when counting, by default 4.

    Returns
    -------
    dict[SmallMoleculeComponent, int|np.int64]
        Dict of nodes (SmallMoleculeComponent) and their corresponding number of cycles counted.
    """
    graph = nx.DiGraph(ligand_network.graph).to_undirected()
    # TODO: check if there is a possibility in nx to omit cycles going back over already chosen nodes (1-2-3-2-1)
    # TODO: if possible remove the corresponding code
    raw_cycles = [
        n
        for c in nx.simple_cycles(graph, length_bound=higher_bound)
        for n in c
        if len(set(c)) == len(c)
    ]
    # TODO: `cou` is currently np.int64 - we should make a true int for consistency
    uni, cou = np.unique(raw_cycles, return_counts=True)

    # Add 0 nodes
    res = dict(zip(uni, cou))
    for n in graph.nodes:
        if n not in res:
            res[n] = 0
    return res


def get_transformation_failure_robustness(
    ligand_network: LigandNetwork,
    failure_rate: float = 0.05,
    nrepeats: int = 100,
    seed: int | None = None,
) -> np.float64:  # TODO: make this a float
    """Estimate the robustness of a LigandNetwork by:
        1. Removing a number of edges corresponding to the percentage given by `failure_rate`.
        2. Randomly repeating step 1 `nrepeats` times.
        3. Average of the results of the repeats in step 2 for a result between 0 (was always disconnected) and 1 (never was disconnected).


    Parameters
    ----------
    ligand_network : LigandNetwork

    failure_rate : float, optional
        Failure rate of edges to use in the estimation, by default 0.05.
    nrepeats : int, optional
        Number of repeat samples to use in the estimation, by default 100.
    seed : int | None, optional
        Seed to use for the random number generator, by default None.

    Returns
    -------
    np.float64
        Average probability that the network remained connected when `failure_rate` fraction of edges was removed.
        Between 0.0 (the network was always disconnected) and 1.0 (network was never disconnected).
    """
    edges = list(ligand_network.edges)
    npics = max(int(np.round(len(edges) * failure_rate)), 1)

    connected = []
    rng = np.random.default_rng(seed=seed)
    for _ in range(nrepeats):
        nn = ligand_network
        edges_to_del = [edges[e] for e in rng.choice(len(edges), npics, replace=False)]
        nn = tools.delete_transformation(nn, edges_to_del, must_stay_connected=False)
        connected.append(float(nn.is_connected()))

    # np.mean on list of bools give expected float answer
    r = np.mean(connected)
    return r
