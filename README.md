<p align="center">
    <picture align="center">
      <source media="(prefers-color-scheme: dark)" srcset="https://github.com/OpenFreeEnergy/konnektor/blob/36fc908f89777b8d67ce837a354adc699de6f405/.img/konnektor_logo_style5.png">
      <source media="(prefers-color-scheme: light)" srcset="https://github.com/OpenFreeEnergy/konnektor/blob/36fc908f89777b8d67ce837a354adc699de6f405/.img/konnektor_logo_style4.png">
      <img alt="Konnektor`s fancy logo" src="https://github.com/OpenFreeEnergy/konnektor/blob/36fc908f89777b8d67ce837a354adc699de6f405/.img/konnektor_logo_style4.png" width=35% >
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

Konnektor offers at the moment basic network planers. 

More will be here soon!
### Implemented Simple Network Layouts
![](.img/network_layouts.png)

### Tools for constructing more complex Networks
Additionally we are working on Molecule diversity based Networks here with HIF2A, which also allows concatenating Networks:
![](.img/Divesity_network.png)


## Usage
```python3
import numpy as np
from openfe_benchmarks import benzenes
from kartograf import KartografAtomMapper
from konnektor.visualization import draw_ligand_network
from openfe.setup.atom_mapping.lomap_scorers import default_lomap_score

compounds = list(filter(lambda x: not x.name in ["lig_2", "lig_3", "lig_4", "lig_7"],
                        benzenes.get_system().ligand_components))

from konnektor.network_planners import (MaximalNetworkPlanner, RadialLigandNetworkPlanner,
                                        MinimalSpanningTreeLigandNetworkPlanner, CyclicLigandNetworkPlanner)

networkers = [MaximalNetworkPlanner, RadialLigandNetworkPlanner,
              MinimalSpanningTreeLigandNetworkPlanner, CyclicLigandNetworkPlanner, DiversityNetworkPlanner]

networks = []
for networker_cls, name in zip(networkers,["Max", "Radial", "MST", "Cyclic", "Div"]):
    networker = networker_cls(mapper=KartografAtomMapper(), scorer=default_lomap_score)
    network = networker.generate_ligand_network(compounds)
    network.name=name
    networks.append(network)

#Visualize
fig, axes = plt.subplots(ncols=2, nrows=3, figsize=[16,3*9])
axes= np.array(axes).flat
fs=22
for ax, net in zip(axes, [max_network, radial_network, mst_network, cyclic_network]):
    draw_ligand_network(network=net, title=net.name, ax=ax, node_size=1500, fontsize=fs)
    ax.axis("off")

axes[-1].axis("off")

fig.show()
```


## Installation
you can install Kartograf via the package manager of your choice:

For Development purposes:
```shell
git clone https://github.com/OpenFreeEnergy/konnektor.git

cd konnektor
conda env create -f environment.yml

conda activate konnektor
pip install .
```
