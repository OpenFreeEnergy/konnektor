================
Network Planning
================

Network planning turns a set of ligands into a concrete plan of which transformations to compute.
**konnektor** does this through two kinds of planner:
**Generators**, which build a network from a set of components, and **Concatenators**,
which join networks that already exist.
In both, each edge is represented by an ``AtomMapping``, which defines the relationship between the two ligands.
An ``AtomMappingScorer`` (any function that takes an ``AtomMapping`` and returns a float in [0,1]) assigns an edge weight proportional to the expected difficulty of that transformation.
The planner then combines those scores with a graph-construction algorithm to choose which edges make up the network.

.. image:: ../_static/img/networks.png

Network Generators
__________________

Network Generators are planners that construct networks from a set of components.
They are usually the starting point for any network planning efforts and come in a wide variety of layouts.

.. image:: ../_static/img/generator.png


konnektor provides Generators across the spectrum from minimal to fully connected:

- **Minimal layouts**: the Star and Minimal Spanning Tree (MST) networks use the fewest edges that still connect the set (N−1),
  so they are cheapest but most sensitive to failures.
- **Redundant layouts**: the Twin Star, Redundant MST, N-Node-Edges and Cyclic networks add extra edges to survive some
  transformation failures, trading cost for robustness.
- **Maximal layouts**: the Maximal Network computes every possible edge, typically as a starting point that is then reduced;
  the Heuristic Maximal Network approximates it with fewer edges.

.. image:: ../_static/img/network_layouts.png


For example, to build a minimal spanning tree network from a set of components:

.. code-block:: python

    from konnektor.utils import toy_data
    from konnektor.network_planners import MinimalSpanningTreeNetworkGenerator

    components, mapper, scorer = toy_data.build_random_dataset(n_compounds=8)

    planner = MinimalSpanningTreeNetworkGenerator(mappers=mapper, scorer=scorer)
    network = planner.generate_ligand_network(components)


Network Concatenators
______________________

Where a Generator builds a network from components, a **Concatenator** joins networks
that already exist.
It applies when two or more networks share no ligands, so there is no common node to
merge on. When the networks *do* share nodes, you can use the merge method, see
:doc:`network_tools`.
The Concatenator instead connects otherwise disjoint networks by introducing new edges between them.

One application is the repair of disconnected ligand networks.
When some transformations in a network fail, the surviving (successful) edges can leave the
ligands split into disconnected subnetworks that can no longer be ranked against one another.
A Concatenator can introduce new edges between these subnetworks to restore connectivity.

.. image:: ../_static/img/concatenator.png


The following example constructs two independent networks and connects them using an
:class:`MstConcatenator`:

.. code-block:: python

    from konnektor.utils import toy_data
    from konnektor.network_planners import MstConcatenator, MinimalSpanningTreeNetworkGenerator

    components, mapper, scorer = toy_data.build_random_dataset(n_compounds=8)
    planner = MinimalSpanningTreeNetworkGenerator(mappers=mapper, scorer=scorer)

    # for illustration: two networks, sharing no ligands
    net_a = planner.generate_ligand_network(components[:4])
    net_b = planner.generate_ligand_network(components[4:])

    concatenator = MstConcatenator(mappers=mapper, scorer=scorer)
    network = concatenator.concatenate_networks([net_a, net_b])   # networks -> one network

For two networks, candidate mappings are generated between ligands belonging to
different networks.
As with a Generator, each candidate ``AtomMapping`` is evaluated using the supplied ``scorer``, and
the Concatenator chooses which connecting edges to keep.

When more than two subnetworks are concatenated, :class:`MstConcatenator` treats each
subnetwork as a node in a graph. A minimum spanning tree over these
subnetworks determines the connections required to produce a connected ligand
network. Consequently, connecting ``k`` subnetworks requires ``k - 1`` new edges.

A single edge connecting two parts of a network is a point of failure.
Where additional redundancy is required, :class:`RedundantMstConcatenator` provides more than one connection per join.

konnektor currently provides three Concatenators:

- :class:`MaxConcatenator`: keeps *every* possible connecting edge. Typically used as a starting point that is then reduced.
- :class:`MstConcatenator`: selects the edges required to connect the subnetworks
  using a minimum spanning tree, giving the minimal reconnection (``k - 1`` edges for ``k`` subnetworks).
- :class:`RedundantMstConcatenator`: constructs multiple (``n_redundancy``) spanning trees, where possible,
  to introduce additional connections and reduce dependence on individual
  connecting edges.

All Concatenators accept ``exclude_edges`` in :func:`concatenate_networks`. These
mappings identify ligand pairs that should be excluded when proposing new
connections. This is useful when repairing a network because transformations that
have already failed can be prevented from being proposed again.

See :doc:`network_tools` for what else you can do with existing networks.
