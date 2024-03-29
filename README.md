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

Konnektor is a package offering multiple ways on constructing free energy networks and modifying those.

👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷

**WARNING: This repository is highly under development and more will be here soon!** 

👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷👷



Try our interactive demo: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/OpenFreeEnergy/konnektor/blob/main/examples/konnektor_example.ipynb#scrollTo=GU32PaMkzD7x)

## Content
### Implemented Simple Network Layouts
Several Network layouts are implemented in Konnektor. Many of them bring 
advantages and disadvantages, check below.
![](.img/Network_Overview.png)

### Tools for constructing more complex Networks
we are implementing tools for more complex Network generation:
![image](https://github.com/OpenFreeEnergy/konnektor/assets/12428005/5fbb253c-f0d3-41bf-bd92-f520b1363b6d)

### Enable More Complex Higher Order Networks
Additionally we are working on Molecule diversity based Networks here with HIF2A, which also allows concatenating Networks:
![image](https://github.com/OpenFreeEnergy/konnektor/assets/12428005/c4ee0b63-7580-4825-b0cb-dc076e4cb9f4)

## Code Example
```python3
# Here we generate some input data.
from openfe_benchmarks import benzenes
compounds = list(filter(lambda x: not x.name in ["lig_2", "lig_3", "lig_4", "lig_7"],
                        benzenes.get_system().ligand_components))

# Pick your Favourite Network layout with favourite AtomMapper and Scorer
from openfe.setup import KartografAtomMapper, lomap_scorers
from konnektor.network_planners import CyclicLigandNetworkPlanner
networker = CyclicLigandNetworkPlanner(mapper=KartografAtomMapper(), 
                                       scorer=lomap_scorers.default_lomap_score)

# Generate Network
network =networker.generate_ligand_network(compounds)
network.name="Cyclic Network"

# Visualize the generated network
from konnektor.visualization import draw_ligand_network
fig = draw_ligand_network(network=network, title=network.name)

fig.show()
```
![example fig](.img/example_out.png)


## Installation

### Latest release
Konnektor can be installed via the package following package managers:

```shell
pip install konnnektor
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
