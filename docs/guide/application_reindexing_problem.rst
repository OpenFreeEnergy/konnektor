==============================================================
Network Plannnig for Free Energy Calculations in Drug Design
==============================================================

In this text we would like to use a classical hands on drug design example use-case, that is searching for the best small molecule to inhibit a protein.
Usually drug design project teams and the computational chemists generate a ton of molecule designs, form which they need to decide on which are relevant and which not.
This can of course be done experimentally by testing different properties of the molecules.
However in order to be more efficient and faster, it is of interest to use computational methods to destinguish between relevant and non relevant molecules.
The computational chemist usually applies a workflow of methods, that start with very fast methods, that sort out the obvious molecules (100 thousands) going to more computational expensive methods, that
are more accurate, but can only be used on smaller amounts of molecules (hundreds).

Relative free energy calculation methods using alchemical MD simulations are more on the more computational expensive end of methods.
Therefore efficient calculation planning is of interest in order to not waste time and compute resources.
But how are can these calculations be made efficient?
This of course is a multi-layered question, starting from methodological aspects like improving the sampling speed or the caclulation time of a timestep by using parallel approaches on GPUs.
But also on a higher level, the required amount of calculations can be adjusted.
This is the level on where Konnektor comes into play and helps you to adjust this..

turorial: :doc:`/tutorial/building_networks`.
NetworkPlanner: :doc:`/guide/network_planner`.
