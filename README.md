[![Logo](https://img.shields.io/badge/OSMF-OpenFreeEnergy-%23002f4a)](https://openfree.energy/)
[![build](https://github.com/OpenFreeEnergy/konnektor/actions/workflows/ci.yaml/badge.svg)](https://github.com/OpenFreeEnergy/konnektor/actions/workflows/ci.yaml)
[![coverage](https://codecov.io/gh/OpenFreeEnergy/konnektor/branch/main/graph/badge.svg)](https://codecov.io/gh/OpenFreeEnergy/konnektor)
[![Documentation Status](https://readthedocs.org/projects/konnektor/badge/?version=latest)](https://konnektor.readthedocs.io/en/latest/?badge=latest)

# Konnektor: The tools for building networks for your FE calcuations
====================================================================

Konnektor offers at the moment basic network planers. 

More will be here soon!

## Usage
```python3
import numpy as np
from openfe_benchmarks import benzenes
from kartograf import KartografAtomMapper
from konnektor.visualization import draw_ligand_network
from openfe.setup.atom_mapping.lomap_scorers import default_lomap_score

#Get Input Data
compounds = list(filter(lambda x: not x.name in ["lig_2", "lig_3", "lig_4", "lig_7"], 
                        benzenes.get_system().ligand_components))

#Build Networks
from konnektor.network_planners import (MaximalNetworkPlanner, RadialLigandNetworkPlanner, 
                                        MinimalSpanningTreeLigandNetworkPlanner, CyclicLigandNetworkPlanner)

networkers = [MaximalNetworkPlanner, RadialLigandNetworkPlanner,
              MinimalSpanningTreeLigandNetworkPlanner, CyclicLigandNetworkPlanner]

networks = []
for networker_cls, name in zip(networkers,["Max", "Radial", "MST", "Cyclic"]):
    networker = networker_cls(mapper=KartografAtomMapper(), scorer=default_lomap_score)
    network = networker.generate_ligand_network(compounds)
    network.name=name
    networks.append(network)

#Visualize
fig, axes = plt.subplots(ncols=2, nrows=2, figsize=[16,9])
axes= np.array(axes).flat
for ax, net in zip(axes, [max_network, radial_network, mst_network, cyclic_network]):
    draw_ligand_network(network=net, title=net.name, ax=ax, node_size=1500)
    ax.axis("off")
    
fig.show()
```

![](.img/network_layouts.png)


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