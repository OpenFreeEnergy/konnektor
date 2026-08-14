=============
Network Tools
=============
Where a planner builds a network, network tools are the
utilities *around* that step. They help you prepare, improve, and maintain free energy
networks, acting on different objects at different points in the workflow:

- :doc:`clustering`: group components by chemical similarity before planning, so you can
  plan within groups rather than across the whole set.
- :doc:`Intermediate generation <intermediate_generators>`: when a single transformation
  is too difficult to compute directly, insert one or more new intermediate ligands
  between two components, splitting one hard edge into easier ones.
- :doc:`network handling <network_handling>`: edit a network you already have.

.. toctree::
   :maxdepth: 1

   clustering
   intermediate_generators
   network_handling
