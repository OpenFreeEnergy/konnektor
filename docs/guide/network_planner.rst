==============================================================
Network Planners
==============================================================

Network Generators
__________________
Building a drug candidate ranking problem is actually a graph/network construction problem,
beacuse a connected network, with the candidates as nodes, can describe all the relations of each candidate to each other.
The relation of two molecules here equals an edge of this network, translating to one RBFE calculation approach.
So in practice computational chemists use many RBFE calculations to build up networks, describing the relations of all drug candidates with each other.

.. image:: ../_static/img/networks.png

Note as the free energy is a thermodynamic state function, it is path independent.
This means that a network only needs to be connected to give insight into all relations (if all RBFE calculations have very good quality).
Examplified, a molecule A and a molecule C have a direct relation x.
But if they are connected via molecule B, using relation y and relation z, the sum of y+z still will yield x.

Ok cool, what can we do with this, as mentioned the FE network should be efficiently calculated.
Keep in mind computational chemists use high performance computing resources in order to compute RBFE calculation, taking several hours.
So the RBFE calculations are expensive and usually it is avoided to calculate to many RBFE calculations.

This can be reflected directly in the FE network design.
In order to orchestrate the calculations, a plan for which connections to calculated is generated prior to the simulations.
This is where Konnektor enters the stage.

Konnektor is a package that supports you to generate the FE network calculation plan.
It implements multiple network layouts, that have different advantages and diadvantages.
How is such a network plan generated?
For our example case with the RBFE calculations, each edge can be realized as an AtomMapping indicating a common substructure between the two molecules, to be compared.
This `AtomMapping`, can then be scored using an `AtomMappingScorer`, indicating how difficult a transformation is expected to be in terms of convergence or accuracy.

An network planning algorithm now can use these scores together with graph constructing algorithms, in order to find the best calculation paths.
The minimal amount of edges for a ranking can be achieved with the Star Network and the MST Network (N-1).

Network Concatenators
______________________
Network Concatenators are solving the problem on how to connect two networks, that do not share any edges.
This can be done by solving a bi-partite graph matching problem.



Network Tools
______________

* Basic Operations

* Analysis

* Visualizations
