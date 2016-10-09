Features
========

Supported Quantities
--------------------

If present in the output of a calculation, the following quantities are recognised by *iago*. All quantities here
consist of raw data as calculated by the simulation package itself. None of them is derived (except for unit
conversions).

=================== ==============
Property            CP2K QuickStep
=================== ==============
MD Step Number      yes
MD Time             yes
Simulation box      yes
Temperature         yes
Pressure            yes
Conserved Quantity  yes
Core Self Energy    yes
Core Hamiltonian    yes
Hartree Energy      yes
XC Energy           yes
HFX Energy          yes
Dispersion Energy   yes
Total Energy        yes
Potential Energy    yes
Kinetic Energy      yes
Temperature Drift   yes
IASD                yes
S2 determinant      yes
SCF Cycles Needed   yes
OT Cycles Needed    yes
ERI Count           yes
=================== ==============

Supported Calculations
----------------------

Based on the trajectory data and the simulation package log file, the following quantities can be calculated.

.. py:currentmodule:: iago.Analyser

.. autosummary::

	Analyser.dynamic_plane
	Analyser.dynamic_distance

