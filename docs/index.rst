Welcome to konnektor's documentation!
=====================================

Konnektor is a Python library for generating networks for free energy calculations.
It contains algorithms and tools for network planning that make setting up the calculation plans much easier.

If you're looking for a tool to perform free energy calculations, check out `openfe <https://docs.openfree.energy/en/latest/>`_, which uses **konnektor**!

Why might you need access to different network generation approaches?
As an example, imagine you are given a set of drug candidates that to be ranked with relative binding free energies.
In theory, you could calculate *all* the possible network transformations to get your ligand ranking (we call this a Maximal Network).
However this leads to an explosion in time and compute cost, so we need more efficient ways on how to calculate a drug candidate ranking.

From a thermodynamic perspective, all the transformations in a Maximal Network are actually required to retrieve a ranking.
In fact, you only need one connection per small molecule to the others, such that you have a connected network, in order to get the ranking.
Examples of this minimal case include the Star Network and Minimal Spanning Tree (MST) Network.
However, these very efficient networks are sensitive to transformation failures, and so network algorithms that add a degree of redundancy are used to improve the network's robustness.
**konnektor** enables you to construct and analyze networks that are custom to your calculations' needs.


In addition to the network planners included in the library, **konnektor** has tooling for:
   - Network modification, such as concatenating networks or deleting edges.
   - Network analysis, including calculating graph scores, getting the connectivities of network nodes, or calculating the network robustness.
   - Network visualization tools and an interactive network visualization widget for use in IPython environments, such as jupyter notebooks.

If you're not ready to install **konnektor** locally, you can explore its functionality with the `interactive demo <https://colab.research.google.com/github/OpenFreeEnergy/konnektor/blob/main/examples/konnektor_example.ipynb#scrollTo=GU32PaMkzD7x>`_ directly in your browser.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   tutorials/index
   guide
   api
   CHANGELOG
