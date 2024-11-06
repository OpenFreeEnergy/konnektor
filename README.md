<p align="center">
    <picture align="center">
      <source media="(prefers-color-scheme: dark)" srcset="https://github.com/OpenFreeEnergy/konnektor/blob/2e52eadb32a4a2813bc81fc3d0499f7563c8b4e6/.img/konnektor_logo_5.png">
      <source media="(prefers-color-scheme: light)" srcset="https://github.com/OpenFreeEnergy/konnektor/blob/2e52eadb32a4a2813bc81fc3d0499f7563c8b4e6/.img/konnektor_logo_4.png">
      <img alt="Konnektor`s fancy logo" src="https://github.com/OpenFreeEnergy/konnektor/blob/2e52eadb32a4a2813bc81fc3d0499f7563c8b4e6/.img/konnektor_logo_4.png" width=35% >
    </picture>
</p>


Konnektor: Tools for Networks in your FE Calculations
====================================================================

[//]: # (Badges)
[![Logo](https://img.shields.io/badge/OSMF-OpenFreeEnergy-%23002f4a)](https://openfree.energy/)
[![build](https://github.com/OpenFreeEnergy/konnektor/actions/workflows/ci.yaml/badge.svg)](https://github.com/OpenFreeEnergy/konnektor/actions/workflows/ci.yaml)
[![coverage](https://codecov.io/gh/OpenFreeEnergy/konnektor/branch/main/graph/badge.svg)](https://codecov.io/gh/OpenFreeEnergy/konnektor)
[![Documentation Status](https://readthedocs.org/projects/konnektor/badge/?version=latest)](https://konnektor.readthedocs.io/en/latest/?badge=latest)

[![Pip Install](https://img.shields.io/badge/pip%20install-konnektor-d9c4b1)](https://pypi.org/project/konnektor/)

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

Checkout our article on Konnektor in the Journal of Chemical Modelling and Information: [𝐊𝐨𝐧𝐧𝐞𝐤𝐭𝐨𝐫: 𝐀 𝐅𝐫𝐚𝐦𝐞𝐰𝐨𝐫𝐤 𝐟𝐨𝐫 𝐔𝐬𝐢𝐧𝐠 𝐆𝐫𝐚𝐩𝐡 𝐓𝐡𝐞𝐨𝐫𝐲 𝐭𝐨 𝐏𝐥𝐚𝐧 𝐍𝐞𝐭𝐰𝐨𝐫𝐤𝐬 𝐟𝐨𝐫 𝐅𝐫𝐞𝐞 𝐄𝐧𝐞𝐫𝐠𝐲 𝐂𝐚𝐥𝐜𝐮𝐥𝐚𝐭𝐢𝐨𝐧𝐬 - Benjamin Ries*, Richard J. Gowers, Hannah M. Baumann, David W. H. Swenson,  Michael M. Henry, James R. B. Eastwood, Irfan Alibay, and David Mobley](https://pubs.acs.org/doi/10.1021/acs.jcim.4c01710).

You can find our pre-print on [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/66d9b574cec5d6c1423643e4).

Try our interactive demo: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/OpenFreeEnergy/konnektor/blob/main/examples/konnektor_example.ipynb#scrollTo=GU32PaMkzD7x)



## Content
### Implemented Network Layouts
Several Network layout generating algorithms are implemented in Konnektor, with different advantages and disadvantages.
From an algorithmic perspective most of the algorithms are actually a reduction method of the Maximal Network.
To speed-up the Maximal Network Algorithm we implemented an parallelization scheme to it. Below you can find some of our layouts:

![image](docs/_static/img/network_layouts.png)

### Tools for handling Networks
Konnektor implements tools, that allow for example to merge (if a node is shared in the networks) or concatenate (if no node is shared) networks,
append single molecules (nodes) to a network or delete transformations/molecules from a network.

![image](docs/_static/img/network_operations_small.png)

### Enable More Complex Higher Order Networks
Another goal of Konnektor is to go beyond the standard network layout algorihtms and allow easy implementation of more complex network algorithms.
This is achieved by combining the Tools and Network Generator Algorithms, to build up to more advanced workflows.

![image](docs/_static/img/advanced_networking.png)

## Code Example

```python3
# Here we generate some input data.
from konnektor.data import get_benzene_ligands

compounds = list(filter(lambda x: not x.name in ["lig_2", "lig_3", "lig_4", "lig_7"],
                        get_benzene_ligands()))

# Pick your Favourite Network layout with favourite AtomMapper and Scorer
from openfe.setup import KartografAtomMapper, lomap_scorers
from konnektor.network_planners import CyclicNetworkGenerator

networker = CyclicNetworkGenerator(mappers=KartografAtomMapper(),
                                   scorer=lomap_scorers.default_lomap_score)

# Generate Network
network = networker.generate_ligand_network(compounds)
network.name = "Cyclic Network"

# Visualize the generated network
from konnektor.visualization import draw_ligand_network
fig = draw_ligand_network(network=network, title=network.name)

fig.show()
```
![example fig](.img/example_out.png)


## Installation

### Latest release
Konnektor can be installed from conda-forge:

```shell
mamba install -c conda-forge konnnektor
```

### Developement version
The developing setup of Konnektor works like this:

```shell
git clone https://github.com/OpenFreeEnergy/konnektor.git

cd konnektor
mamba env create -f environment.yml

mamba activate konnektor
pip install -e .

```

## License
This library is made available under the MIT open source license.

## Authors

The OpenFE development team.

## Acknowledgments
Thanks to Enrico Ruijsenaars, Jenke Scheen and Josh Horton for great discussions!
