Welcome to konnektor's documentation!
=====================================

Konnektor is a Python library for generating networks for free energy calculations.
It contains algorithms and tools for network planning that make setting up the calculation plans much easier.

If you're looking for a tool to perform free energy calculations, check out `openfe <https://docs.openfree.energy/en/latest/>`_, which uses **konnektor**!.

Why might you need access to different network generation approaches?
As an example, imagine you are given a set of drug candidates that to be ranked with relative binding free energies.
In theory, you could calculate _all_ the possible network transformations to get your ligand ranking (we call this a Maximal Network).
However this leads to an explosion in time and compute cost, so we need more efficient ways on how to calculate a drug candidate ranking.

From a thermodynamic perspective not all transformations are actually required to retrieve a ranking.
In fact, you only need one connection per small molecules to the others in order to get the ranking, for example a Star Network or Minimal Spanning Tree (MST) Network.
However, these very efficient networks are sensitive to transformation failures, and so network algorithms that add a degree of redundancy improves the calculation's robustness.

In addition to the network planners included in the library, **konnektor** has tooling for:
   - Network modification, such as concatenating networks or deleting edges.
   - Network analysis, including calculating graph scores, getting the connectivities of network nodes, or calculating the network robustness.
   - Network visualization tools and an interactive network visualization widget for use in IPython environments, such as jupyter notebooks.

Try our interactive demo: |Colab|

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   tutorials/index
   guide
   api
   CHANGELOG

.. |Colab| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/OpenFreeEnergy/konnektor/blob/main/examples/konnektor_example.ipynb#scrollTo=GU32PaMkzD7x
