=====================
Installation
=====================

User Setup
=============

Kartograf will be packaged with OpenFE but we need more time for that. :)

Alternatively, you can install Kartograf also as a standalone package via pip
or soon conda:

``pip install konnektor``



Developer Setup
================

For Developers, we recommend to setup a Konnektor environment like in the
following example, allowing modification of the code live in the package:

    git clone https://github.com/OpenFreeEnergy/konnektor.git

    cd konnektor
    mamba env create -f environment.yml

    mamba activate konnektor
    python -m pip install -e .

Happy coding! :)
