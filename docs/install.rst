============
Installation
============

**konnektor** comes installed with **openfe** as of ``openfe v1.5.0``!

If you want to use **konnektor** as a standalone package, see below.

Installation from conda-forge
=============================

**konnektor** can be installed from conda-forge using your choice of micromamba (recommended), mamba, or conda:

```shell
micromamba create -c -conda-forge -n konnektor konnektor
micromamba activate konnektor

```

Developer Installation
======================

For developers, we recommend setting up a local developer installation so that your changes to the code are immediately reflected in the functionality.


First, clone the git repo:

```shell
    git clone https://github.com/OpenFreeEnergy/konnektor.git
    cd konnektor
```

Then create and activate a conda environment that includes all of **konnektor's** dependencies:

```shell
    mamba env create -f environment.yaml
    mamba activate konnektor
```

Finall, create an editable installation:

```shell
    python -m pip install -e .
```

Happy coding!
