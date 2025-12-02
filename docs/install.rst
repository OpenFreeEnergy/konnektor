=====================
Installation
=====================

User Setup
=============

Konnektor will be packaged with OpenFE, soon. :)

Konnektor can be installed via the package following package managers:

```shell
mamba -c -conda-forge install konnektor
```

Developer Setup
================

For Developers, we recommend to setup a Konnektor environment like in the
following example, allowing modification of the code live in the package:

    git clone https://github.com/OpenFreeEnergy/konnektor.git

    cd konnektor
    mamba env create -f environment.yaml

    mamba activate konnektor
    python -m pip install -e .

Happy coding! :)
