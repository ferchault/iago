This document describes the database structure used for any DatabaseProvider subclass. In any case, there is exactly one database for each bucket which several tables as outlined below

=========== ================== =======================================
Table name  Human readable     Typical contents
=========== ================== =======================================
meta        meta data          Information that is valid for the whole bucket, i.e. comments. Has exactly two columns: *key*, *value*.
config      input files        Implemented as a ordered tree structure where each node can have a scalar or vector value associated. First tree nodes are the run names.
run         per run            Information that is valid for a single run, i.e. timing information.
frame       per frame          Information that is valid for a single frame, i.e. set of distances, energies, box dimensions. 
atom        per atom           Information that is always valid for a single atom for all runs, i.e. atom types, elements.
=========== ================== =======================================

Information is stored at the hierarchical level at which it is produced. E.g. an explicitly configured ensemble will show up in the *config* tree as well as in the *run* table while the implicit default value for an ensemble will only be visible in the *run* table.
