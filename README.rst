Terminology
===========

================= ====================
Name              Meaning
================= ====================
Run               A single molecular dynamics (MD) run or single point calculation.
Bucket            One set of runs that belong together. In particular, they share the physical system.
Frame             Single configuration of the system. May be part of a time series (MD) or not (single point).
================= ====================


Data Storage Levels
===================

In atomistic simulations, data can be related to different entities. iago deals with the following:

- Per bucket and run (e.g. configuration values)
- Per atom and bucket (e.g. atom kind, element, atom name)
- Per atom and frame (e.g. coordinates, velocities, forces)
- Per frame and bucket (e.g. total energy)
