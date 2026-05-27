Why use **konnektor** ?
=======================

Why might you need access to different network generation approaches?
As an example, imagine you are given a set of drug candidates that to be ranked with relative binding free energies.
In theory, you could calculate *all* the possible network transformations to get your ligand ranking (we call this a Maximal Network).
Though robust, a Maximal Network approach leads to explosion in time and compute cost, and so more efficient networks are needed.

From a thermodynamic perspective, not all the transformations in a Maximal Network are actually required to retrieve a ranking.
In fact, the opposite extreme - a minimally connected network such as a Star Network or a Minimal Spanning Tree (MST) Networks - is actually needed to compute rankings.
However, these very efficient networks are highly sensitive to transformation failures, and so network algorithms that add a degree of redundancy are needed to improve the network's robustness.
**konnektor** enables you to construct and analyze the multitude of possible networks that fall between these extremes to find an appropriate network generation scheme for a given set of ligands.

See the next section for how to get started generating ligand networks with **konnektor**.
