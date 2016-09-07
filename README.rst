===========
Terminology
===========

================= ====================
Name              Meaning
================= ====================
Run               A single molecular dynamics (MD) run or single point calculation.
Bucket            One set of runs that belong together. In particular, they share the physical system.
Frame             Single configuration of the system. May be part of a time series (MD) or not (single point).
================= ====================
===================
Data Storage Levels
===================

In atomistic simulations, data can be related to different entities. iago deals with the following:

- Per bucket and run (e.g. configuration values)
- Per atom and bucket (e.g. atom kind, element, atom name)
- Per atom and frame (e.g. coordinates, velocities, forces)
- Per frame and bucket (e.g. total energy)
============
Config Files
============

There are two locations that hold configuration files: the analysis machine where a reference for the different locations is needed and the location remote where raw, aggregated and metadata is available. In any case, the file format follows the Windows INI files as implemented by the python module `ConfigParser <https://docs.python.org/2/library/configparser.html>`_.

----------------
Analysis Machine
----------------
This file holds the configuration options where to look for buckets on both local and remote machines. An example illustrates the supported protocols.

::

  [localmachine]
  url=file:///home/username/data/
  [remoteserver]
  url=ssh://username@hostname:/path
  skip=True

This file must reside in the user's home directory under the name *.iago.conf*.
