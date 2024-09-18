==============================================================
Network Plannnig for Free Energy Calculations in Drug Design
==============================================================

In this section, we will explore a classical hands-on drug design example,
specifically the search for the best small molecule to inhibit a protein for drug development.
Typically, drug design project teams and computational chemists generate a large
number of molecule designs from which they must determine relevance.

While experimental testing can evaluate various properties of these molecules,
it is can be more efficient to use computational methods to distinguish between
relevant and irrelevant candidates.
The computational chemist usually follows a workflow that starts with very fast
methods to filter out the obvious molecules (in the range of hundreds of thousands),
then progresses to more computationally expensive methods that are more accurate
but can only be applied to smaller sets of molecules (usually in the hundreds).

Relative free energy calculation methods using alchemical MD simulations fall on
the more computationally expensive side of the spectrum. Therefore, making the calculations
efficient is crucial to avoid wasting time and computing resources.
But how can these calculations be made efficient?

This is a multi-layered question, starting with methodological aspects such as
improving sampling speed or reducing timestep calculation time through parallel
approaches on GPUs. At a higher level, the required number of calculations can
also be adjusted. This is where Konnektor comes into play, helping you optimize this process.

In general, there is a trade-off to consider in the calculation network
between efficiency regarding the number of edges and redundancy within the
network. Ideally, a minimal number of edges, such as in a Star Network or
Minimal Spanning Tree (MST) Network, would make calculations most efficient.
However, sometimes an edge calculation may fail due to poor edge selection,
leading to a disconnected network, which results in an incomplete molecule
ranking.

These failures can arise from factors that were not accounted for during the
selection process, such as unforeseen issues or computational outages, like
a node crash. To address this, redundancy can help ensure immediate
resolution of disconnections. Additionally, redundancy can be utilized to
estimate simulation errors through cycle closure analysis, which can also
enhance the final results.

More redundant networks include the Twin Star Network, Redundant MST
Network, N-Node Edges Network, and Cyclic Network. Notably, the Twin Star
Network and the Cyclic Network focus on constructing a set of cycles.

So in applied cases we suggest the slightly more redundant networks in order to build up a robust network for the calculations.

turorial: :doc:`/tutorial/building_networks`.
NetworkPlanner: :doc:`/guide/network_planner`.
