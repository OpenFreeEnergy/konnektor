Welcome to Konnektors's documentation!
=========================================

Konnektor is a package supporting you in planning your free calculations.
It contains multiple algorithms and tools for network planning, that make setting up the calculation plans much easier.
As an example imagen you are given a set of drug candidates that shall be ranked with relative binding free energies.
In theory you could calculate all the possible network transformations, in order to get your ligand ranking (we call this a Maximal Network).
However this leads to an explosion in time and compute cost, therefore we need more efficient ways on how to caluclate a drug candidate ranking.
From a thermodynamic perspective not all transformations are actually required to retrieve a ranking.
In fact you only need one conection per small molecules to the others in order to get the ranking, like for example in Star Networks or Minimal Spanning Tree (MST) Networks.
However we found the very efficient networks to be sensitive to transformation failures, this can be solved with network building algorithms, that are slightly more redundant.

Ontop of the described ligand network planners, Konnektor gives access to tools, that allow for example to concatenate networks or delete transformations of a network.
Analysis of networks, like calculating graph scores, getting the connectivities of network nodes or calculating the network robustness are available too.
Last we want to bring to your attention our Network visualization tools and the provided interactive network visualization widget for IPython like in Jupyter-Lab/Notebooks.

Try our interactive demo: |Colab|

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   tutorial
   guide
   api
   CHANGELOG

.. |Colab| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/OpenFreeEnergy/konnektor/blob/main/examples/konnektor_example.ipynb#scrollTo=GU32PaMkzD7x
