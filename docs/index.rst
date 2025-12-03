Welcome to konnektor's documentation!
=====================================

**konnektor** is a Python library for planning, modifying, and analyzing free energy transformation networks.

If you're looking for a tool to perform free energy calculations, check out `openfe <https://docs.openfree.energy/en/latest/>`_, which uses **konnektor**!

Why might you need access to different network generation approaches?
As an example, imagine you are given a set of drug candidates that to be ranked with relative binding free energies.
In theory, you could calculate *all* the possible network transformations to get your ligand ranking (we call this a Maximal Network).
Though robust, a Maximal Network approach leads to explosion in time and compute cost, and so more efficient networks are needed.

From a thermodynamic perspective, not all the transformations in a Maximal Network are actually required to retrieve a ranking.
In fact, the opposite extreme - a minimally connected network such as a Star Network or a Minimal Spanning Tree (MST) Networks - is actually needed to compute rankings.
However, these very efficient networks are highly sensitive to transformation failures, and so network algorithms that add a degree of redundancy are needed to improve the network's robustness.
**konnektor** enables you to construct and analyze the multitude of possible networks that fall between these extremes to find an appropriate network generation scheme for a given set of ligands.

In addition to network planning algorithms, **konnektor** includes tooling for:
   - Network modification, such as concatenating networks or deleting edges.
   - Network analysis, including calculating graph scores, getting the connectivities of network nodes, or calculating the network robustness.
   - Network visualization tools and an interactive network visualization widget for use in IPython environments, such as jupyter notebooks.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   tutorials/index
   guide
   api
   CHANGELOG
