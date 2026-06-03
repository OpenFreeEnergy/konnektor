===============
Node Clustering
===============
The `Clusterer` classes can be used to separate a collection of Components by a specific property into sub-groups.
They serve as powerful tools within a workflow and enable some network generation schemes, such as the Starry Sky Network.

In **konnektor**, the following clustering classes are currently available:

* `ChargeClusterer`: separates molecules by net charge changes.
* `ScaffoldClusterer`: separates molecules by shared scaffolds
* `DiversityClusterer`: uses fingerprints to cluster the different molecules.
