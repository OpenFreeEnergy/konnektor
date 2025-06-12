=========
Changelog
=========

.. current developments

v0.2.0
====================

**Changed:**

* For both ``ExplicitNetworkGenerator.generate_network_from_indices()`` and ``ExplicitNetworkGenerator.generate_network_from_names()``, all ``components`` will be included as nodes in the network. If this results in a disconnected network (in the case where one or more of the ``components`` is not included in ``indices``, or ``names``, respectively), a warning will be raised (`#122 <https://github.com/OpenFreeEnergy/konnektor/pull/122>`_).
* If multiple mappers are passed to a network planner, but no scorer, the *first* mapping from the *first* mapper provided will be used. This is a change from the previous behavior, which would use the *first* mapping from the *last* mapper. The current behavior is now consistent with ``openfe``\'s network planning behavior (`#84 <https://github.com/OpenFreeEnergy/konnektor/pull/84>`_).
* ``clustering.auxilliary_featurizer.ChargeTransformer`` parameter ``parallel`` renamed to ``n_jobs`` which allows for designating the maximum number of concurrently running jobs (`#133 <https://github.com/OpenFreeEnergy/konnektor/pull/133>`_).
* For radial (star) network generation, if ``central_component`` is already in the list passed to ``components``, a self-edge will no longer be created, and a warning is raised (`#145 <https://github.com/OpenFreeEnergy/konnektor/pull/145>`_).
* **konnektor** now treats mapping scores of 1.0 as best, and 0.0 as worst. This is an inversion from the prior behavior (0.0 best, 1.0 worst) to match the behavior of the rest of the openfe ecosystem (`#138 <https://github.com/OpenFreeEnergy/konnektor/pull/138>`_).
* Renamed ``RadialLigandNetworkPlanner`` to ``RadialNetworkGenerator`` to be consistent with the other network generator class names (`#160 <https://github.com/OpenFreeEnergy/konnektor/pull/160>`_).
* For ``RedundantMinimalSpanningTreeNetworkGenerator``, if the number of redundancies (``n_redundancy``) is larger than the number of redundant MSTs able to be generated for the network, a warning will be raised, and the resulting network will contain the most redundant MSTs able to be generated (`#122 <https://github.com/OpenFreeEnergy/konnektor/pull/122>`_).
